#!/usr/bin/python3
"""A python module to simplify supervision and running of multiple DynamO simulations.

You define state variables which you'd like to vary during a parameter
sweep as well as output properties you'd like to collect. This module
then manages running the simulations, sorting them if the state
variables have changed/increased/etc, continuing the simulations if
they were aborted or the parameters were changed, as well as
calculating error estimates for the properties obtained at each state
point.
"""

import bz2
import glob
import math
# Include everything "standard" in here. Try to keep external
# dependencies only imported when they are used, so this can be
# easily deployed on a cluster.
import os
import subprocess
import sys
import time
from multiprocessing import Pool, cpu_count

import alive_progress
import scipy

#import xml.etree.cElementTree as ET

try:
    import lxml.etree as ET
    print("Running with lxml.etree")
except ImportError:
    try:
        # normal cElementTree install
        import xml.etree.cElementTree as ET
        print("Running with cElementTree")
    except ImportError:
        try:
            # normal ElementTree install
            import xml.etree.ElementTree as ET

            import elementtree.ElementTree as ET
            print("Running with ElementTree")
        except ImportError:
            print("Failed to import ElementTree from any known place")

import numpy as np
import uncertainties
from ./datastat import WeightedArray, WeightedFloat, linear_interp


class SkipThisPoint(BaseException):
    pass

def print_to_14sf(f):
    """Utility function to print a variable to 14 significant figures"""
    if isinstance(f, str):
        return f
    return '{:.{p}g}'.format(f, p=14)

def conv_to_14sf(f):
    if isinstance(f, float):
        return float(print_to_14sf(f))
    else:
        return f

class XMLFile:
    """A wrapper around  to allow loading and saving to
    bzip2 compressed files. """
    
    def __init__(self, filename, compressed=None):
        """Loads the xml file, decompressing it with bz2 first if the
        filename ends with .bz2"""
        import io
        if isinstance(filename, str):
            self._filename = filename        
            if  filename.endswith('.xml.bz2'):
                import bz2
                f = bz2.BZ2File(filename)
                self.tree = ET.parse(f)
                f.close()
            elif filename.endswith('.xml'):
                self.tree = ET.parse(filename)
            else:
                raise RuntimeError('Unknown file extension for configuration file load "'+filename+'"')
        elif isinstance(filename, io.BytesIO):
            data = filename.read()
            try:
                self._filename = filename.name
            except:
                self._filename = "bytestream.xml"
            try:
                import bz2
                data = bz2.decompress(data)
                print("Successfully decompressed!")
            except:
                print("Could not decompress, trying direct reading")
            self.tree = ET.fromstring(data)
        else:
            raise RuntimeError("Could not determine file type")
        
    def save(self, filename):
        if filename.endswith('.xml.bz2'):
            import bz2
            f = bz2.BZ2File(filename, mode='w')
            f.write(ET.tostring(self.tree.getroot()))
            f.close()
        elif filename.endswith('.xml'):
            open(filename, 'w').write(ET.tostring(self.tree.getroot()))
        else:
            raise RuntimeError('Unknown file extension for configuration file save "'+filename+'"')

    def __str__(self):
        return "XMLFile("+self._filename+")"

# A XMLFile/ElementTree but specialised for DynamO configuration files
class ConfigFile(XMLFile):
    def __init__(self, filename):
        super().__init__(filename)

    # Number of particles in the config
    def N(self):
        return int(len(self.tree.findall('.//Pt')))

    # Primary image volume
    def V(self):
        dims = self.image_dimensions()
        return dims[0] * dims[1] * dims[2]

    # number density
    def n(self):
        return self.N() / self.V()
    
    def __getitem__(self, key):
        return ConfigFile.config_props[key]['recalc'](self)
    
    def __str__(self):
        return "ConfigFile("+self._filename+")"

    def histogramTether(self, limits=[None, None, None]):
        tethers = self.tree.findall('.//Global/CellOrigins/Origin')
        particles = self.tree.findall('.//Pt/P')
        data = np.ndarray((len(tethers), 3))
        
        for idx, tether, particle in zip(range(len(tethers)), tethers, particles):
            p = np.array(list(map(float, [particle.attrib['x'], particle.attrib['y'], particle.attrib['z']])))
            t = np.array(list(map(float, [tether.attrib['x'], tether.attrib['y'], tether.attrib['z']])))
            data[idx] = p-t
        return np.histogramdd(data, range=[limits, limits, limits], bins=11)

    def histogramTether1D(self, limits=[None]):
        tethers = self.tree.findall('.//Global/CellOrigins/Origin')
        particles = self.tree.findall('.//Pt/P')
        data = np.ndarray((len(tethers), 1))
        for idx, tether, particle in zip(range(len(tethers)), tethers, particles):
            p = np.array(list(map(float, [particle.attrib['x'], particle.attrib['y'], particle.attrib['z']])))
            t = np.array(list(map(float, [tether.attrib['x'], tether.attrib['y'], tether.attrib['z']])))
            data[idx] = math.sqrt((p-t).dot(p-t))
        return np.histogramdd(data, bins=11)

    def image_dimensions(self):
        V = self.tree.find('.//SimulationSize')
        return [float(V.attrib['x']), float(V.attrib['y']), float(V.attrib['z'])]
    
    def to_freud(self):
        import freud
        N = self.N()
        frame = np.zeros((N, 3), dtype=np.float32)
        box = self.image_dimensions()
        box = freud.box.Box(Lx = box[0], Ly = box[1], Lz = box[2])
        particles = self.tree.findall('.//Pt/P')
        
        for idx, particle in enumerate(particles):
            frame[idx, 0] = np.float32(particle.attrib['x'])
            frame[idx, 1] = np.float32(particle.attrib['y'])
            frame[idx, 2] = np.float32(particle.attrib['z'])
        
        return box, frame
    
    config_props = {}

# A XMLFile/ElementTree but specialised for DynamO output files
class OutputFile(XMLFile):
    def __init__(self, filename):
        super().__init__(filename)
        
    def N(self):
        return int(self.tree.find('.//ParticleCount').attrib['val'])

    def events(self):
        return int(self.tree.find('.//Duration').attrib['Events'])

    def t(self):
        return float(self.tree.find('.//Duration').attrib['Time'])

    def numdensity(self):
        return float(self.tree.find('.//Density').attrib['val'])

    def __getitem__(self, key):
        return OutputFile.output_props[key](self)

    def __str__(self):
        return "OutputFile("+self._filename+")"
    
    output_props = {}

def validate_xmlfile(filename):
    try:
        ET.parse(bz2.BZ2File(filename))
        return True
    except Exception as e:
        print("#!!!#", filename, e)
        return False
    
def validate_outputfile(filename):
    return validate_xmlfile(filename)
    
def validate_configfile(filename):
    return validate_xmlfile(filename)
    
import pickle as pickle


#This function actually sets up and runs the simulations and is run in parallel
def worker(state, workdir, outputplugins, particle_equil_events, particle_run_events, particle_run_events_block_size, setup_worker):
    try:
        if True:
            if not os.path.isdir(workdir):
                #Make the work directory
                os.mkdir(workdir)
                #Save the state
                pickle.dump(state, open(os.path.join(workdir, "state.pkl"), 'wb'))
        
            logfile = open(os.path.join(workdir, 'run.log'), 'a')
            
            print("\n", file=logfile)
            print("################################", file=logfile)
            print("#        Setup Config          #", file=logfile)
            print("################################  ", file=logfile, flush=True)
        
            startconfig = os.path.join(workdir, "start.config.xml.bz2")
            if not os.path.isfile(startconfig) or not validate_configfile(startconfig):
                print("No (valid) config found, creating...", file=logfile, flush=True)
                try:
                    setup_worker(startconfig, state, logfile, particle_equil_events)
                except SkipThisPoint as e:
                    #Leave the work dir, we'll just skip the point
                    return
                except subprocess.CalledProcessError as e:
                    raise RuntimeError('Failed while running setup worker, command was\n"'+str(e.cmd)+'"\nSee logfile "'+str(os.path.join(workdir, 'run.log'))+'"')
            else:
                print("Initial config found.", file=logfile, flush=True)
        
            #Do the equilibration run
            inputfile = startconfig
            outputfile = os.path.join(workdir, '0.config.xml.bz2')
            datafile = os.path.join(workdir, '0.data.xml.bz2')
        
            #Parse how many particles there are
            inconfig = ConfigFile(inputfile)
            N = inconfig.N()
            inconfig = None # Just doing this to free the XML reader
            
            print("\n", file=logfile)
            print("################################", file=logfile)
            print("#      Equilibration Run       #", file=logfile)
            print("################################\n", file=logfile, flush=True)
            
            from subprocess import check_call

            #Only actually do the equilibration if the output data/config is missing
            if not os.path.isfile(outputfile) or not validate_configfile(outputfile) or not os.path.isfile(datafile) or not validate_outputfile(datafile):
                check_call(["dynarun", inputfile, '-o', outputfile, '-c', str(N * particle_equil_events), "--out-data-file", datafile], stdout=logfile, stderr=logfile)
            else:
                print("Found existing valid equilibration run", file=logfile)
        
            #Now do the production runs
            counter = 1
            curr_particle_events = 0
            while curr_particle_events < particle_run_events:
                print("\n", file=logfile)
                print("################################", file=logfile)
                print("#        Production Run        #", file=logfile)
                print("################################", file=logfile, flush=True)
                print("Events ",curr_particle_events, "/", particle_run_events, "\n", file=logfile, flush=True)
                inputfile = os.path.join(workdir, str(counter-1)+'.config.xml.bz2')
                # Abort if input file is missing
                if not os.path.isfile(inputfile):
                    print("ERROR! input file missing?", file=logfile)
                    return
                outputfile = os.path.join(workdir, str(counter)+'.config.xml.bz2')
                datafile = os.path.join(workdir, str(counter)+'.data.xml.bz2')
                dotherun = False
                if (not os.path.isfile(outputfile)):
                    print("output config file for run "+str(counter)+" "+outputfile+" is missing, doing the run", file=logfile)
                    dotherun = True
                elif (not validate_configfile(outputfile)):
                    print("output config file for run "+str(counter)+" "+outputfile+" is corrupted, doing the run", file=logfile)
                    dotherun = True
                elif (not os.path.isfile(datafile)):
                    print("output data file for run "+str(counter)+" "+datafile+" is missing, doing the run", file=logfile)
                    dotherun = True
                elif (not validate_outputfile(datafile)):
                    print("output data file for run "+str(counter)+" "+datafile+" is corrupted, doing the run", file=logfile)
                    dotherun = True

                if dotherun:
                    check_call(["dynarun", inputfile, '-o', outputfile, '-c', str(N * particle_run_events_block_size), "--out-data-file", datafile]+outputplugins, stdout=logfile, stderr=logfile)
                    curr_particle_events += particle_run_events_block_size
                    counter += 1
                else:
                    of = OutputFile(datafile)
                    events_per_N_run = of.events() / of.N()
                    curr_particle_events += events_per_N_run
                    print("Found existing config and data for run "+str(counter)+" with "+str(events_per_N_run)+"N events, skipping", file=logfile)
                    counter += 1
        
                #Process the output data now
            print("\n", file=logfile)
            print("################################", file=logfile)
            print("#        Run Complete          #", file=logfile)
            print("################################", file=logfile)
            print("Events ",curr_particle_events, "/", particle_run_events, "\n", file=logfile, flush=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError('Failed while running worker, command was\n"'+str(e.cmd)+'"\nSee logfile "'+str(os.path.join(workdir, 'run.log'))+'"')
        
def perdir(args):
    output_dir, particle_equil_events, manager = args
    output_dir = os.path.join(manager.workdir, output_dir)
    if not os.path.isdir(output_dir):
        return {}
    
    configs = glob.glob(os.path.join(output_dir, "*.config.xml.bz2"))
    if len(configs) == 0:
        return {}

    try:
        state = pickle.load(open(os.path.join(output_dir, "state.pkl"), 'rb'))
    except FileNotFoundError:
        print("Skipping: Could not find state.pkl in ", output_dir)
        return {}
    statedict, state = make_state(state)

    #Filter to only the set states (if enabled)
    if manager.only_current_statevars and state not in manager.states:
        return {}
        
    counter = 0
    executed_events = 0

    dataout = {}
    while True:
        if True:
        #try:
            configfilename = os.path.join(output_dir, str(counter)+'.config.xml.bz2')
            datafilename = os.path.join(output_dir, str(counter)+'.data.xml.bz2')
            counter += 1
            
            if (not os.path.isfile(configfilename)) or (not os.path.isfile(datafilename)):
                break
            
            outputfile = OutputFile(datafilename)
            run_events = outputfile.events()
            N = outputfile.N()
            
            if executed_events < particle_equil_events * N:
                executed_events += run_events
                continue
            
            if "NEventsTot" not in dataout:
                dataout["NEventsTot"] = 0
            dataout["NEventsTot"] += outputfile.events()
                
            if "tTotal" not in dataout:
                dataout["tTotal"] = 0
            dataout["tTotal"] += outputfile.t()
                            
            for prop in manager.outputs:
                outputplugin = OutputFile.output_props[prop]
                result = outputplugin.result(state, outputfile, configfilename, counter, manager, output_dir)
                if result != None:
                    if prop not in dataout:
                        dataout[prop] = outputplugin.init()
                    dataout[prop] += result
        #except Exception as e:
        #    print("Processing", output_dir, " gave exception", e)
        #    #raise
    return {state: dataout}

def make_state(state):
    #Convert anything (list, tuple, dict) to a state dictionary
    statedict = dict(state)
    #Then make a tuple of tuples (so its hashable) and sort so the order is indifferent to user choice
    return statedict, tuple(sorted([(statevar, statedict[statevar]) for statevar in statedict]))    
    
import shutil


def reorg_dir_worker(args):
    entry, manager = args
    oldpath = os.path.join(manager.workdir, entry)
    if os.path.isdir(oldpath):
        configs = glob.glob(os.path.join(oldpath, "*.config.xml.bz2"))
        if len(configs) == 0:
            return []
        oldstate = pickle.load(open(os.path.join(oldpath, "state.pkl"), 'rb'))
        oldstatedict, oldstate = make_state(oldstate)

        XMLconfig = ConfigFile(configs[0])
        newstate = oldstatedict.copy()
        for statevar in manager.used_statevariables:
            can_regen = ConfigFile.config_props[statevar]['recalculable']
            regen_val = XMLconfig[statevar]
            # If we can regenerate the state from the
            # config file, then do that to verify the
            # state value. If not, then only replace the
            # state value if it is missing from the old
            # state.
            if can_regen or (statevar not in newstate):
                newstate[statevar] = regen_val

        newstatedict, newstate = make_state(newstate)
        if  newstate != oldstate:
            #print("Would rewrite oldstate ", repr(oldstate), "to", repr(newstate))
            pickle.dump(newstate, open(os.path.join(oldpath, "state.pkl"), 'wb'))
        
        return [(oldpath, newstate)]
    return []

class SimManager:
    def __init__(self, workdir, statevars, outputs, restarts=1, processes=None):
        if not shutil.which("dynamod"):
            raise RuntimeError("Could not find dynamod executable.")

        if not shutil.which("dynarun"):
            raise RuntimeError("Could not find dynamod executable.")

        self.restarts = restarts
        self.outputs = set(outputs)
        self.workdir = workdir
        self.output_plugins = set()

        if len(statevars) == 0:
            raise RuntimeError("We need some state variables to work on")
        
        #Make sure the state vars are in ascending order for easy output to screen
        self.statevars = [[(key, sorted(value)) for key, value in sweep] for sweep in statevars]
       
        #Now create all states for iteration. We need to do this now,
        #as we need to get the generated state variables too.
        self.states = self.iterate_state(self.statevars)

        #We need a list of the state variables used 
        self.used_statevariables = list(map(lambda x : x[0], next(iter(self.states))))
        
        for output in self.outputs:
            # Check output is defined
            if output not in OutputFile.output_props:
               raise RuntimeError('The "'+output+'" output is not defined!') 
            
            # Check the outputs have all the required statevars
            for dep_statevar in OutputFile.output_props[output]._dep_statevars:
                if dep_statevar not in self.used_statevariables:
                    raise RuntimeError('The "'+output+'" output requires the "'+dep_statevar+'" state variable, but its missing from state vars')
            
            # Add any dependent outputs
            for dep_output in OutputFile.output_props[output]._dep_outputs:
                self.outputs.add(dep_output)
        
            #Collect what output plugins are needed
            for dep_outputplugin in OutputFile.output_props[output]._dep_outputplugins:
                self.output_plugins.add(dep_outputplugin)

        self.output_plugins = list(self.output_plugins)
            
        if not os.path.isdir(workdir):
            os.mkdir(workdir)
        self.processes = processes
        if self.processes is None:
            self.processes = cpu_count()

    def getstatedir(self, state, idx):
        return os.path.join(self.workdir, self.statename(state) + "_" + str(idx))
    
    def statename(self, state, var_separator='_'):
        output = ""
        #When building the name, we make sure state vars are ordered.
        statedict = dict(state)
        for statevar, stateval in [(statevar, statedict[statevar]) for statevar in self.used_statevariables]:
            if isinstance(stateval,float):
                stateval = print_to_14sf(stateval)
            else:
                stateval = str(stateval)
            output = output + statevar + "_" + stateval + var_separator
        return output[:-1]
    
    def getnextstatedir(self, state, oldpath = None):
        idx = 0
        while True:
            newpath = self.getstatedir(state, idx)
            if oldpath == newpath:
                return newpath
            if not os.path.isdir(newpath) and not os.path.isfile(newpath):
                return newpath            
            idx += 1
            
    def reorg_dirs(self):
        entries = os.listdir(self.workdir)
        print("Reorganising existing data directories...")
        n = len(entries)
        
        pool = Pool(processes=self.processes)
        with alive_progress.alive_bar(n) as progress:
            #This is a parallel loop, returning items as they finish in arbitrary order
            for actions in pool.imap_unordered(reorg_dir_worker, [(d, self) for d in entries], chunksize=10):
                #We do the actual moving here. We only want one thread
                #doing the moving as it must check what directories
                #already exist to figure out the new name. This is
                #hard to do in parallel
                for oldpath, newstate in actions:
                    newpath = self.getnextstatedir(newstate, oldpath=oldpath)
                    if oldpath != newpath:
                        shutil.move(oldpath, newpath)
                progress()
        
    def iterate_state(self, statevars):
        # Loop over all permutations of the state variables. We use a
        # set to automatically remove repeats, especially as we're
        # running these states in parallel.
        states = set()
        
        import itertools
        for sweep in statevars:
            sweep_states = 1
            #Make a list of state variables, and a list of their values
            statevar, statevals = zip(*sweep)

            #Then combine all permutations of the state values 
            for stateval in itertools.product(*statevals):
                #Use the values to build something that can be used as a dictionary
                state = {var: conv_to_14sf(val) for var, val in zip(statevar, stateval)}
                for svar,sval in list(state.items()): #We slice to make sure modifying the original doesn't derail the loop
                    if 'gen_state' in ConfigFile.config_props[svar]:
                        state = ConfigFile.config_props[svar]['gen_state'](state)
                _, state = make_state(state)
                states.add(state)
        return states

    def get_run_files(self, workdir, min_events, max_events=None):
        if max_events is None:
            max_events = float("inf")
        counter = 0
        curr_particle_events = 0

        #First, get past the min_events configs
        equil_configs=[]
        while curr_particle_events < min_events:
            outputfile = os.path.join(workdir, str(counter)+'.config.xml.bz2')
            datafile = os.path.join(workdir, str(counter)+'.data.xml.bz2')

            #Check both files exist, if not, bail!
            if not os.path.isfile(outputfile) or not validate_configfile(outputfile) or not os.path.isfile(datafile) or not validate_outputfile(datafile):
                return equil_configs, []
            equil_configs.append((outputfile, datafile))
            of = OutputFile(datafile)
            counter += 1
            curr_particle_events += of.events() / of.N()

        run_configs=[]
        while curr_particle_events < max_events:
            outputfile = os.path.join(workdir, str(counter)+'.config.xml.bz2')
            datafile = os.path.join(workdir, str(counter)+'.data.xml.bz2')

            #Check both files exist, if not, bail!
            if not os.path.isfile(outputfile) or not validate_configfile(outputfile) or not os.path.isfile(datafile) or not validate_outputfile(datafile):
                return equil_configs, run_configs
            
            run_configs.append((outputfile, datafile))
            of = OutputFile(datafile)
            counter += 1
            curr_particle_events += of.events() / of.N()
            
                        
    def run(self, setup_worker, particle_equil_events, particle_run_events, particle_run_events_block_size):            
        print("Generating simulation tasks for the following sweeps")
        for idx, sweep in enumerate(self.statevars):
            print(" Sweep", idx)
            for statevar, statevals in sweep:
                print("  ",statevar, "∈", list(map(print_to_14sf, statevals)))

        tot_states = len(self.states)

        running_tasks = []
        tasks_completed = 0
        tasks_failed = 0
        task_count = 0
        errors = []
        pool = Pool(processes=self.processes)
        
        print("Building task tree...")
        class Task:
            def __init__(self, workertuple):
                self._workertuple = workertuple
                self._next_tasks = []

            def is_done(self):
                return self._result.ready()

            def is_successful(self):
                return self._result.successful()
            
            def start(self):
                self._result = pool.apply_async(worker, args=self._workertuple)
                
            def to_follow(self, task):
                self._next_tasks.append(task)

            def next_tasks(self):
                return self._next_tasks

            def failed(self):
                return 1 + sum([task.failed() for task in self._next_tasks])
        
        #We break up tasks into blocks of events
        for state in self.states:
            for idx in range(self.restarts):
                workdir = self.getstatedir(state, idx)
                run_events = 0
                parent_task = None
                while run_events < particle_run_events:
                    run_events += particle_run_events_block_size
                    new_task = Task((state, workdir, self.output_plugins, particle_equil_events, run_events, particle_run_events_block_size, setup_worker))
                    if parent_task is None:
                        running_tasks.append(new_task)
                    else:
                        parent_task.to_follow(new_task)
                    parent_task = new_task
                    task_count += 1

        print("Running", len(self.states) * self.restarts, "state points as ", task_count, "simulation tasks in parallel with", self.processes, "processes")

        with alive_progress.alive_bar(task_count, manual=True) as progress:
            for task in running_tasks:
                task.start()

            ok = True
            while len(running_tasks) > 0:
                still_running = []
                for task in running_tasks:
                    if not task.is_done():
                        still_running.append(task)
                    else:
                        if task.is_successful():
                            tasks_completed += 1
                            if ok:
                                for nxttask in task.next_tasks():
                                    nxttask.start()
                                    still_running.append(nxttask)
                        else:
                            tasks_completed += task.failed()
                            ok = False
                            try:
                                task._result.get()
                            except Exception as e:
                                print("\n ERROR: Found error ", e, "\n")
                                errors.append(e)
                
                running_tasks = still_running
                progress(tasks_completed / task_count)
                time.sleep(0.5)

        print("Terminating and joining threads...")
        pool.close()
        pool.join()
        
        if len(errors) > 0:
            import traceback
            print("\nERROR: Found",len(errors),"exceptions while processing")
            print(''.join(traceback.format_exception(errors[0])))
            f=open("error.log", 'w')
            print(''.join([''.join(traceback.format_exception(error)) for error in errors]), file=f)
            
            print('Remaining errors written to "error.log"')
            raise RuntimeError("Parallel execution failed")

    def fetch_data(self, particle_equil_events, only_current_statevars = False):
        self.only_current_statevars = only_current_statevars
        
        output_dirs = os.listdir(self.workdir)
        print("Fetching data...")
        n = len(output_dirs)

        pool = Pool(processes=self.processes)

        import collections

        #We store the extracted data in a dict of dicts. The first
        #dict is for the state, the second for the property.
        state_data = collections.defaultdict(dict)        

        #So we run the per data dir operation, then reduce everything
        state_data = {}
        with alive_progress.alive_bar(n) as progress:
            #This is a parallel loop, returning items as they finish in arbitrary order
            for result in pool.imap_unordered(perdir, [(d, particle_equil_events, self) for d in output_dirs], chunksize=10):
                #Here we process the returned data from a single directory
                for state, data in result.items():
                    if state not in state_data:
                        state_data[state] = data
                    else:
                        target = state_data[state]
                        for key, value in data.items():
                            if key not in target:
                                target[key] = value
                            else:
                                target[key] += value
                progress()
        
        #We now prep the data for processing, we convert our
        #WeightedFloat's to ufloats as pandas supports that natively.
        for statevars, data in state_data.items():
            for prop in data:
                if isinstance(data[prop], WeightedFloat):
                    data[prop] = data[prop].ufloat()

        #Here we create the dataframe
        import pandas
        df = pandas.DataFrame([{**dict(state), **output} for state, output in state_data.items()])

        
        #If there's no data, then there's no columns so the next
        #bit fails. Avoid that
        if len(df) == 0:
            return df
        
        #Here, we're just adjusting the column order to follow what was given by the user.
        cols = list(df.columns.values)
        for statevar in self.used_statevariables:
            cols.remove(statevar)
        df = df[[statevar for statevar in self.used_statevariables]+cols]
        #Now we sort items by the state variables in the order given.
        df = df.sort_values(by=[statevar for statevar in self.used_statevariables])

        #Now we write out the data
        import pickle
        pickle.dump(df, open(self.workdir+".pkl", 'wb'))

        return df

# ###############################################
# #          Definition of state vars           #
# ###############################################
ConfigFile.config_props["N"] = {'recalculable':True, 'recalc': lambda config: config.N()}
ConfigFile.config_props["ndensity"] = {'recalculable':True, 'recalc': lambda config: conv_to_14sf(config.n())}
ConfigFile.config_props["InitState"] = {'recalculable':False, 'recalc': lambda config: "FCC"}


def Lambda_config(XMLconfig):
    tag = XMLconfig.tree.find('.//Interaction[@Type="SquareWell"]')
    if tag is None:
        return float('inf')
    else:
        return conv_to_14sf(float(tag.attrib['Lambda']))
ConfigFile.config_props["Lambda"] = {'recalculable':True, 'recalc': lambda config: Lambda_config}

def Rso_config(XMLconfig):
    tag = XMLconfig.tree.find('.//Global[@Type="SOCells"]')
    if tag is None:
        return float('inf')
    else:
        return conv_to_14sf(float(tag.attrib['Diameter']) / 2)    
ConfigFile.config_props["Rso"] = {'recalculable':True, 'recalc':Rso_config}


def PhiT_config(XMLconfig):
    Rso = Rso_config(XMLconfig)
    if Rso == float('inf'):
        return float('inf')
    density = XMLconfig.n()
    phiT = density * math.pi * (4.0/3.0) * (Rso**3)
    return conv_to_14sf(phiT)
def PhiT_gen(state):
    density = state["ndensity"]
    phiT = state["PhiT"]
    if phiT == float('inf'):
        Rso = float('inf')
        if 'Rso' in state:
            if state['Rso'] != float('inf'):
                raise RuntimeError("State contains conflicting Rso and PhiT")
        else:
            state['Rso'] = float('inf')
    else:
        Rso = (phiT / (density * math.pi *(4.0/3.0))) ** (1/3.0)
        if 'Rso' in state:
            if state['Rso'] != conv_to_14sf(Rso):
                raise RuntimeError("State contains conflicting Rso and PhiT")
        else:
            state['Rso'] = conv_to_14sf(Rso)
    return state
ConfigFile.config_props["PhiT"] = {'recalculable':True, 'recalc':PhiT_config, 'gen_state':PhiT_gen}

def kT_config(XMLconfig):
    tag = XMLconfig.tree.find('.//System[@Name="Thermostat"]')
    if tag is None:
        return float('inf')
    else:
        return conv_to_14sf(float(tag.attrib['Temperature']))
ConfigFile.config_props["kT"] = {'recalculable':True, 'recalc':kT_config}


class OutputProperty:
    def __init__(self, dependent_statevars : list, dependent_outputs : list, dependent_outputplugins: list):
        self._dep_statevars = dependent_statevars
        self._dep_outputs = dependent_outputs
        self._dep_outputplugins = dependent_outputplugins

    def init(self):
        return None

    def result(self, state, outputfile, configfilename, counter, manager, output_dir):
        return None
    
class SingleAttrib(OutputProperty):
    def __init__(self, tag, attrib, dependent_statevars, dependent_outputs, dependent_outputplugins, time_weighted=True, div_by_N=False, div_by_t=False, missing_val = 0, skip_missing=False):
        OutputProperty.__init__(self, dependent_statevars, dependent_outputs, dependent_outputplugins)
        self._tag = tag
        self._attrib = attrib
        self._time_weighted = time_weighted
        self._div_by_N = div_by_N
        self._div_by_t = div_by_t
        self._missing_val = missing_val
        self._skip_missing=skip_missing

    def init(self):
        return WeightedFloat()

    def value(self, outputfile):
        tag = outputfile.tree.find('.//'+self._tag)
        if tag is None:
            if self._missing_val is None:
                raise RuntimeError('Failed to find the tag "'+self._tag+'" in the outputfile')
            else:
                return self._missing_val
            
        if self._attrib not in tag.attrib:
            if self._missing_val is None:
                raise RuntimeError('Failed to find attribute "'+self._attrib+'" in the tag "'+self._tag+'" in the outputfile')
            else:
                return self._missing_val
                
        val = float(tag.attrib[self._attrib])
        if self._div_by_N:
            N = outputfile.N()
            if N > 0:
                val /= outputfile.N()
        if self._div_by_t:
            t = outputfile.t()
            if t > 0:
                val /= t
        return val

    def weight(self, outputfile):
        tag = outputfile.tree.find('.//'+self._tag)
        #If we have a missing_val defined, then we can use and weight it, else don't give this any weight
        if ((self._missing_val is None) or self._skip_missing) and (tag is None or self._attrib not in tag.attrib):
            return 0
        
        if self._time_weighted:
            return float(outputfile.tree.find('.//Duration').attrib['Time'])
        else:
            return float(outputfile.tree.find('.//Duration').attrib['Events'])

    def result(self, state, outputfile, configfilename, counter, manager, output_dir):
        return WeightedFloat(self.value(outputfile), self.weight(outputfile))

def parseToArray(text):
    data = []
    for row in text.split('\n'):
        row_data = list(map(float, row.split()))
        if len(row_data) > 0:
            data.append(row_data)
    return np.array(data)

class CollisionMatrixOutputProperty(OutputProperty):
    def __init__(self):
        OutputProperty.__init__(self, dependent_statevars=[], dependent_outputs=[], dependent_outputplugins=['-LCollisionMatrix'])

    def result(self, state, outputfile, configfilename, counter, manager, output_dir):
        return None

class VACFOutputProperty(OutputProperty):
    def __init__(self):
        OutputProperty.__init__(self, dependent_statevars=[], dependent_outputs=[], dependent_outputplugins=['-LVACF'])

    def result(self, state, outputfile, configfilename, counter, manager, output_dir):
        restart_idx = output_dir.split('_')[-1]
        filename_root = manager.workdir+'_VACF/' + manager.statename(state, var_separator='/') + "/run_" + restart_idx + '_' + str(counter)
        os.makedirs(filename_root, exist_ok=True)
        for tag in outputfile.tree.findall('.//VACF/Particles/Species'):
            pickle.dump(parseToArray(tag.text), open(filename_root + '/species_'+tag.attrib['Name']+'.pkl', 'wb'))
        
        for tag in outputfile.tree.findall('.//VACF/Topology/Structure'):
            pickle.dump(parseToArray(tag.text), open(filename_root + '/topology_'+tag.attrib['Name']+'.pkl', 'wb'))
        return None

class RadialDistributionOutputProperty(OutputProperty):
    def __init__(self):
        OutputProperty.__init__(self, dependent_statevars=[], dependent_outputs=[], dependent_outputplugins=['-LRadialDistribution'])

    def init(self):
        return WeightedArray()

    def result(self, state, outputfile, configfilename, counter, manager, output_dir):
        #Presume that each tag is in order, and has a common bin width
        samples = float(outputfile.tree.find('.//RadialDistributionMoments').attrib["SampleCount"])
        bin_width = 0.01        
        #Grab all the moments, drop the R values for now
        moment_tags = outputfile.tree.findall('.//RadialDistributionMoments/Species/Moment')
        moments = np.array([[float(line.split()[1]) for line in tag.text.strip().split("\n")] for tag in moment_tags])

        # In the simulation we have collected the moments of, N(r), the number of pairs below a radius r. 
        # These moments are collected about an origin/offset N₀(r), i.e. ⟨(N(r)-N₀(r))^n⟩
        # The origin is a numerical trick to reduce the size of the moments to reduce precision issues, its an
        # approximation of the mean, by taking the initial value of N(r) at t=0.

        #moments is a 6xNbins array with N0, then <(N-N0)>, <(N-N0)^2>, and <(N-N0)^3> etc.
        N0 = np.copy(moments[0,:])
        #We set the N0 column (copied above) into the <(N-N0)^0> column. This makes the formula exact 
        moments[0,:] = np.ones_like(moments[0,:])

        # We eventually want to convert to the central moments ⟨(N(r)-⟨N(r)⟩)^n⟩, To
        # do this we need to determine the first moment ⟨N(r)⟩. There is a general
        # expression for transforming moments:
        #
        # ⟨(x-b)^n⟩ = Σ_{i=0}^n comb(n, i) ⟨(x-a)^i⟩ (a-b)^{n-i}
        #
        # comb is combination function in scipy.special (binomial coefficient).
        # But the answer for this first step is simple, x=N(r), b=0, a=N₀(r), n=1 gives the straightforward identity
        #
        # ⟨N(r)⟩ = ⟨N(r)-N₀(r)⟩ + N₀(r)

        # For all other steps n>2, b=⟨N(r)⟩, and a = N₀(r).
        # ⟨(N(r)-⟨N(r)⟩)^n⟩ = Σ_{i=0}^n comb(n, i) ⟨(N(r)-N₀(r))^i⟩ (N₀(r)-⟨N(r)⟩)^{n-i}

        central_moments = np.zeros((moments.shape[0]-1, moments.shape[1]))
        central_moments[0,:] = moments[1,:] + N0
        ##Now make the other central moments
        for n in range(2, moments.shape[0]):
            central_moments[n-1] = sum(scipy.special.comb(n, i) * moments[i] * (N0-central_moments[0,:])**(n-i) for i in range(n+1))

        return WeightedArray(central_moments, samples)

class RadialDistEndOutputProperty(OutputProperty):
    def __init__(self):
        OutputProperty.__init__(self, dependent_statevars=[], dependent_outputs=[], dependent_outputplugins=[])


    def result(self, state, outputfile, configfilename, counter, manager, output_dir):
        #This ending is not used, we only want one output for speed.
        #restart_idx = output_dir.split('_')[-1]
        #+ "/run_" + restart_idx + '_' + str(counter)
        filename_root = manager.workdir+'_RadialDist/' + manager.statename(state, var_separator='/')

        if os.path.isdir(filename_root):
            return None

        os.makedirs(filename_root, exist_ok=True)
        
        logfile = open(os.path.join(filename_root, 'run.log'), 'a')
        from subprocess import check_call
        check_call(["dynarun", configfilename, '-o', os.path.join(filename_root, 'RadDist.config.xml.bz2'), '-c', '0', "--out-data-file", os.path.join(filename_root, 'RadDist.out.xml.bz2'), '-LRadialDistribution'], stdout=logfile, stderr=logfile)
        
        of = OutputFile(os.path.join(filename_root, 'RadDist.out.xml.bz2'))
        for tag in of.tree.findall('.//RadialDistribution/Species'):
            A = tag.attrib['Name1']
            B = tag.attrib['Name2']
            output_pkl = filename_root + '/species_'+A+'_'+B+'.pkl'
            pickle.dump(parseToArray(tag.text), open(output_pkl, 'wb'))
        return None

class OrderParameterProperty(OutputProperty):
    '''
    See here https://freud.readthedocs.io/en/stable/modules/order.html#freud.order.Steinhardt
    
    This property is expensive to run at data collection time, as it processes every configuration file to determine the order parameter.
    The advantage is that it can be run on any simulation, no need for extra output plugins.
    '''
    def __init__(self, L):
        OutputProperty.__init__(self, dependent_statevars=[], dependent_outputs=[], dependent_outputplugins=[])
        self.L = L
    
    def result(self, state, outputfile, configfilename, counter, manager, output_dir):
        # We use freud to calculate the Steinhardt order parameter
        # This function is called on every "run" of the simulation, with the outputfile and config file
        # So we process the particle positions in the config file to calculate the order parameter
        import freud
        configfile = ConfigFile(configfilename)
        box, points = configfile.to_freud()

        #Steinhardt for FCC
        ql = freud.order.Steinhardt(self.L)
        ql.compute((box, points), neighbors={"num_neighbors": self.L})
        ql_value = ql.particle_order
        
        return WeightedFloat(np.mean(ql_value), 1)

    def init(self):
        return WeightedFloat()
    
    def result(self, state, outputfile, configfilename, counter, manager, output_dir):
        import freud
        configfile = ConfigFile(configfilename)
        box, points = configfile.to_freud()

        #Steinhardt for FCC
        ql = freud.order.Steinhardt(self.L)
        ql.compute((box, points), neighbors={"num_neighbors": self.L})
        ql_value = ql.particle_order
        
        return WeightedFloat(np.mean(ql_value), 1)

class ChungLuConfigurationModel(OutputProperty):
    '''
    https://en.wikipedia.org/wiki/Configuration_model#Chung-Lu_Configuration_Model
    

    '''
    def __init__(self, L):
        OutputProperty.__init__(self, dependent_statevars=[], dependent_outputs=[], dependent_outputplugins=[])
        self.L = L
    
    def result(self, state, outputfile, configfilename, counter, manager, output_dir):
        
        configfile = ConfigFile(configfilename)
        if len(configfile.tree.findall('.//Interaction/CaptureMap')) != 1:
            raise RuntimeError("ChungLuConfigurationModel requires a single CaptureMap in the config file")

        import networkx as nx
        G = nx.Graph()
        for particle in configfile.tree.findall('.//Pt'):
            G.add_node(int(particle.attrib['ID']))

        for pair in configfile.tree.findall('.//Interaction/CaptureMap/Pair'):
            G.add_edge(int(pair.attrib['ID1']), int(pair.attrib['ID2']))
        
        return WeightedFloat(nx.community.modularity(G), 1)

    def init(self):
        return WeightedFloat()
    
    def result(self, state, outputfile, configfilename, counter, manager, output_dir):
        import freud
        configfile = ConfigFile(configfilename)
        box, points = configfile.to_freud()

        #Steinhardt for FCC
        ql = freud.order.Steinhardt(self.L)
        ql.compute((box, points), neighbors={"num_neighbors": self.L})
        ql_value = ql.particle_order
        
        return WeightedFloat(np.mean(ql_value), 1)


OutputFile.output_props["N"] = SingleAttrib('ParticleCount', 'val', [], [], [], missing_val=None)#We use missing_val=None to cause an error if the tag is missing
OutputFile.output_props["p"] = SingleAttrib('Pressure', 'Avg', [], [], [], missing_val=None)
OutputFile.output_props["cv"] = SingleAttrib('ResidualHeatCapacity', 'Value', [], [], [], div_by_N=True, missing_val=None)
OutputFile.output_props["u"] = SingleAttrib('UConfigurational', 'Mean', [], [], [], div_by_N=True, missing_val=None)
OutputFile.output_props["T"] = SingleAttrib('Temperature', 'Mean', [], [], [], missing_val=None)
OutputFile.output_props["density"] = SingleAttrib('Density', 'val', [], [], [], missing_val=None)
OutputFile.output_props["MSD"] = SingleAttrib('MSD/Species', 'diffusionCoeff', [], [], ['-LMSD'], missing_val=None, skip_missing=True)
OutputFile.output_props["NeventsSO"] = SingleAttrib('EventCounters/Entry[@Name="SOCells"]', # Outputfile tag name
                                                    'Count', # Outputfile tag attribute name
                                                    ["Rso"], # Required state variable
                                                    [], # Required output variables
                                                    [], # Required output plugins
                                                    div_by_N=True, # Divide the count by N
                                                    div_by_t=True, # Also divide by t
                                                    missing_val=0) # If counter is missing, return 0
OutputFile.output_props["VACF"] = VACFOutputProperty()
OutputFile.output_props["RadialDistEnd"] = RadialDistEndOutputProperty()
OutputFile.output_props["RadialDistribution"] = RadialDistributionOutputProperty()
OutputFile.output_props["FCCOrder"] = OrderParameterProperty(6)
OutputFile.output_props["CollisionMatrix"] = CollisionMatrixOutputProperty()

if __name__ == "__main__":

    pass
