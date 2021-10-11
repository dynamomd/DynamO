#!/usr/bin/python3
"""A python module to simplify supervision and running of multiple DynamO simulations.

You define state variables which you'd like to vary during a parameter
sweep as well as output properties you'd like to collect. This module
then manages running the simulations, sorting them if the state
variables have changed/increased/etc, continuing the simulations if
they were aborted or the parameters were changed, as well as
calculating error estimates for the properties obtained at each state
point.

This file uses python's multiprocessing package, so you cannot use
this library in interactive mode unless you pass processes=1 to the
SimManager.run() method.

"""

# Include everything "standard" in here. Try to keep external
# dependencies only imported when they are used, so this can be
# easilly deployed on a cluster.
import os, glob, progress_bars, sys, time, math, subprocess, bz2

from multiprocessing import Pool, TimeoutError, cpu_count

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

import uncertainties
import numpy as np
from datastat import WeightedFloat, linear_interp

class SkipThisPoint(BaseException):
    pass

def print_to_14sf(f):
    """Utility function to print a variable to 14 significant figures"""
    if isinstance(f, str):
        return f
    return '{:g}'.format(float('{:.{p}g}'.format(f, p=14)))


def conv_to_14sf(f):
    if isinstance(f, float):
        return float(print_to_14sf(f))
    else:
        return f

class XMLFile:
    """A wrapper around  to allow loading and saving to
    bzip2 compressed files. """
    
    def __init__(self, filename):
        """Loads the xml file, decompressing it with bz2 first if the
        filename ends with .bz2"""
        self._filename = filename
        
        if filename.endswith('.xml.bz2'):
            import bz2
            f = bz2.BZ2File(filename)
            self.tree = ET.parse(f)
            f.close()
        elif filename.endswith('.xml'):
            self.tree = ET.parse(filename)
        else:
            raise RuntimeError('Unknown file extension for configuration file load "'+filename+'"')
        
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
        V = self.tree.find('.//SimulationSize')
        return float(V.attrib['x']) * float(V.attrib['y']) * float(V.attrib['z'])

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

        
    config_props = {}

def statename(statevars):
    output = ""
    for statevar, stateval in statevars:
        if isinstance(stateval,float):
            stateval = print_to_14sf(stateval)
        else:
            stateval = str(stateval)
        output = output + statevar + "_" + stateval + "_"
    return output[:-1]

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
def worker(state, workdir, idx, outputplugins, particle_equil_events, particle_run_events, particle_run_events_block_size, setup_worker):
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
        
    
import shutil

class SimManager:
    def __init__(self, workdir, statevars, outputs, restarts=1, processes=None):
        self.workdir = workdir
        self.statevars = statevars
        self.statevars_dict = dict(statevars)
        self.restarts = restarts
        self.outputs = set(outputs)
        
        if not shutil.which("dynamod"):
            raise RuntimeError("Could not find dynamod executable.")

        if not shutil.which("dynarun"):
            raise RuntimeError("Could not find dynamod executable.")

        self.output_plugins = set()

        #Make a sample state, using just the first set of variables,
        #to figure out what set AND generated state variables.
        _, test_state = self.iterate_state([(statevar, vals[0:1]) for statevar, vals in self.statevars])
        test_state = test_state[0][0]
        
        for output in self.outputs:
            # Check output is defined
            if output not in OutputFile.output_props:
               raise RuntimeError('The "'+output+'" output is not defined!') 
            
            # Check the outputs have all the required statevars
            for dep_statevar in OutputFile.output_props[output]._dep_statevars:
                if dep_statevar not in test_state:
                    raise RuntimeError('The "'+output+'" output requires the "'+dep_statevar+'" state variable, but its missing!')

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
        return os.path.join(self.workdir, statename(state) + "_" + str(idx))
        
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
        with progress_bars.ProgressBar(sys.stdout) as progress:
            n = len(entries)
            for i,entry in enumerate(entries):
                progress.update(i / n)
                oldpath = os.path.join(self.workdir, entry)
                if os.path.isdir(oldpath):
                    configs = glob.glob(os.path.join(oldpath, "*.config.xml.bz2"))
                    if len(configs) == 0:
                        continue                
                    oldstate = pickle.load(open(os.path.join(oldpath, "state.pkl"), 'rb'))
                    XMLconfig = ConfigFile(configs[0])
                    newstate = oldstate
                    for statevar,staterange in self.statevars:
                        can_regen = ConfigFile.config_props[statevar]['recalculable']
                        regen_val = XMLconfig[statevar]
                        # If we can regenerate the state from the
                        # config file, then do that to verify the
                        # state value. If not, then only replace the
                        # state value if it is missing from the old
                        # state.
                        if can_regen:
                            newstate[statevar] = regen_val
                        elif statevar not in newstate:
                            newstate[statevar] = regen_val
                        
                    newstatelist = [(statevar, newstate[statevar]) for statevar,staterange in self.statevars]
                    newpath = self.getnextstatedir(newstatelist, oldpath=oldpath)

                    if oldpath != newpath:
                        shutil.move(oldpath, newpath)
                    if  newstate != oldstate:
                        pickle.dump(newstate, open(os.path.join(newpath, "state.pkl"), 'wb'))


    def iterate_state(self, statevars):
        # Loop over all permutations of the state variables
        statevar, statevals = zip(*statevars)
        import itertools
        states = itertools.product(*statevals)

        # As itertools.product doesn't support len() we need to figure out the total number of states ourselves
        tot_states = self.restarts
        for s in statevals:
            tot_states = tot_states * len(s)
            
        retval = []
        for i, stateval in enumerate(itertools.product(*statevals)):
            state = [(var, conv_to_14sf(val)) for var, val in zip(statevar, stateval)]

            #For each state variable, calculate any derived state as needed
            newstate = state[:]
            
            for svar,sval in state:
                if 'gen_state' in ConfigFile.config_props[svar]:
                    newstate = ConfigFile.config_props[svar]['gen_state'](newstate)
            
            state = newstate
            for idx in range(self.restarts):
                retval.append((dict(state), self.getstatedir(state, idx), idx))
        return tot_states, retval

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
        print("Generating simulation tasks for the following ranges")
        for statevar, statevals in self.statevars:
            print(" ",statevar, "∈", list(map(print_to_14sf, statevals)))

        tot_states, states = self.iterate_state(self.statevars)

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
        for state, workdir, idx in states:
            #First we
            run_events = 0
            parent_task = None
            while run_events < particle_run_events:
                run_events += particle_run_events_block_size
                new_task = Task((state, workdir, idx, self.output_plugins, particle_equil_events, run_events, particle_run_events_block_size, setup_worker))
                if parent_task is None:
                    running_tasks.append(new_task)
                else:
                    parent_task.to_follow(new_task)
                parent_task = new_task
                task_count += 1

        print("Running", tot_states, "state points as ", task_count, "simulation tasks in parallel with", self.processes, "processes")
        
        with progress_bars.ProgressBar(sys.stdout) as progress:
            progress.update(0)

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
                progress.update(tasks_completed / task_count)
                time.sleep(0.5)

        print("Terminating and joining threads...")
        pool.close()
        pool.join()
        
        if len(errors) > 0:
            import traceback
            print("\nERROR: Found",len(errors),"exceptions while processing")
            print(''.join(traceback.format_exception(etype=type(errors[0]), value=errors[0], tb=errors[0].__traceback__)))
            f=open("error.log", 'w')
            print(''.join([''.join(traceback.format_exception(etype=type(error), value=error, tb=error.__traceback__)) for error in errors]), file=f)
            
            print('Remaining errors written to "error.log"')
            raise RuntimeError("Parallel execution failed")

    def fetch_data(self, particle_equil_events, only_current_statevars=False):
        output_dirs = os.listdir(self.workdir)

        import collections
        #We store the extracted data in a dict of dicts. The first
        #dict is for the state, the second for the property.
        state_data = collections.defaultdict(dict)

        #As we want to look up by the state, we need a hashable type
        #holding the state. This means we need to fix an order in
        #which the state vars will be used as a key. So we just grab an order once and reuse it

        print("Fetching data...")
        with progress_bars.ProgressBar(sys.stdout) as progress:
            n = len(output_dirs)
            for i, output_dir in enumerate(output_dirs):
                progress.update(i / n)
                output_dir = os.path.join(self.workdir, output_dir)
                if os.path.isdir(output_dir):
                    configs = glob.glob(os.path.join(output_dir, "*.config.xml.bz2"))
                    if len(configs) == 0:
                        continue
                    statedict = pickle.load(open(os.path.join(output_dir, "state.pkl"), 'rb'))
                
                    #XMLconfig = ConfigFile(configs[0])
                    state = tuple((statevar, statedict[statevar]) for statevar,staterange in self.statevars)
                
                    # Check the config is part of the running statevars
                    if only_current_statevars:
                        for statevar, val in statedict.items():
                            if val not in self.statevars_dict[statevar]:
                                continue
                
                    dataout = state_data[state]
                    if "NEventsTot" not in dataout:
                        dataout["NEventsTot"] = 0
                    if "tTotal" not in dataout:
                        dataout["tTotal"] = 0

                    counter = 0
                    executed_events = 0
                    while True:
                        try:
                            outputfilename = os.path.join(output_dir, str(counter)+'.config.xml.bz2')
                            datafilename = os.path.join(output_dir, str(counter)+'.data.xml.bz2')
                            counter += 1
                            
                            if (not os.path.isfile(outputfilename)) or (not os.path.isfile(datafilename)):
                                break
                            
                            #if (not validate_configfile(outputfilename)) or (not validate_outputfile(datafilename)):
                            #    break
                            
                            outputfile = OutputFile(datafilename)
                            run_events = outputfile.events()
                            N = outputfile.N()
                            
                            if executed_events < particle_equil_events * N:
                                #If we weren't passed the equil_events
                                #before starting this sim, then
                                #discard the run as equilibration
                                executed_events += run_events
                                continue
                            
                            #This counts just the "production" events and time
                            dataout["NEventsTot"] += outputfile.events()
                            dataout["tTotal"] += outputfile.t()
                            
                            for prop in self.outputs:
                                if prop not in dataout:
                                    dataout[prop] = OutputFile.output_props[prop].init()
                                dataout[prop] += OutputFile.output_props[prop].result(outputfile)
                            #print(state[0][1], OutputFile.output_props[prop].value(outputfile))
                        except Exception as e:
                            print("Processing", output_dir, " gave exception", e)
                            break

            for statevars, data in state_data.items():
                for prop in data:
                    if isinstance(data[prop], WeightedFloat):
                        data[prop] = data[prop].ufloat()
                
        import pandas
        df = pandas.DataFrame([{**dict(state), **output} for state, output in state_data.items()])

        # Don't continue processing if there was no data to process
        if df.empty:
            return df
        
        cols = list(df.columns.values)
        for statevar, staterange in self.statevars:
            cols.remove(statevar)
        df = df[[statevar for statevar, staterange in self.statevars]+cols]
        df = df.sort_values(by=[statevar for statevar, staterange in self.statevars])
        return df

# ###############################################
# #          Definition of state vars           #
# ###############################################
ConfigFile.config_props["N"] = {'recalculable':True, 'recalc': lambda config: config.N()}
ConfigFile.config_props["ndensity"] = {'recalculable':True, 'recalc': lambda config: conv_to_14sf(config.n())}
ConfigFile.config_props["InitState"] = {'recalculable':False, 'recalc': lambda config: "FCC"}
ConfigFile.config_props["Lambda"] = {'recalculable':False, 'recalc': lambda config: 1.5}
#
def Rso_config(XMLconfig):
    tag = XMLconfig.tree.find('.//Global[@Type="SOCells"]')
    if tag is None:
        return float('inf')
    else:
        return conv_to_14sf(float(tag.attrib['Diameter']) / 2)    
ConfigFile.config_props["Rso"] = {'recalculable':True, 'recalc':Rso_config}
#
def PhiT_config(XMLconfig):
    Rso = Rso_config(XMLconfig)
    if Rso == float('inf'):
        return float('inf')
    density = XMLconfig.n()
    phiT = density * math.pi * (4.0/3.0) * (Rso**3)
    return conv_to_14sf(phiT)
def PhiT_gen(state):
    dictstate = dict(state)
    density = dictstate["ndensity"]
    phiT = dictstate["PhiT"]
    Rso = (phiT / (density * math.pi *(4.0/3.0))) ** (1/3.0)
    if 'Rso' in dictstate:
        if dictstate['Rso'] != conv_to_14sf(Rso):
            raise RuntimeError("State contains conflicting Rso and PhiT")
    else:
        state.append(('Rso', conv_to_14sf(Rso)))
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

    def result(self, outputfile):
        return WeightedFloat(self.value(outputfile), self.weight(outputfile))

def parseToArray(text):
    data = []
    for row in text.split('\n'):
        row_data = list(map(float, row.split()))
        if len(row_data) > 0:
            data.append(row_data)
    return np.array(data)
    
class VACFOutputProperty(OutputProperty):
    def __init__(self):
        OutputProperty.__init__(self, dependent_statevars=[], dependent_outputs=[], dependent_outputplugins=['-LVACF'])

    def init(self):
        return []
        
    def result(self, outputfile):
        result = {}
        result['weight'] = int(outputfile.tree.find('.//VACF').attrib['ticks'])
        for tag in outputfile.tree.findall('.//VACF/Particles/Species'):
            result['Species:'+tag.attrib['Name']] = parseToArray(tag.text)
        for tag in outputfile.tree.findall('.//VACF/Topology/Structure'):
            result['Topology:'+tag.attrib['Name']] = parseToArray(tag.text)
        return [result]

class TransportProperty(OutputProperty):
    def __init__(self):
        OutputProperty.__init__(self, dependent_statevars=[], dependent_outputs=[], dependent_outputplugins=['-LMisc'])

    def init(self):
        return []
        
    def result(self, outputfile):
        result = {}
        for tag in outputfile.tree.findall('.//VACF/Particles/Species'):
            result['Species:'+tag.attrib['Name']] = parseToArray(tag.text)
        for tag in outputfile.tree.findall('.//VACF/Topology/Structure'):
            result['Topology:'+tag.attrib['Name']] = parseToArray(tag.text)
        return [result]
    

OutputFile.output_props["N"] = SingleAttrib('ParticleCount', 'val', [], [], [], missing_val=None)#We use missing_val=None to cause an error if the tag is missing
OutputFile.output_props["p"] = SingleAttrib('Pressure', 'Avg', [], [], [], missing_val=None)
OutputFile.output_props["T"] = SingleAttrib('Temperature', 'Mean', [], [], [], missing_val=None)
OutputFile.output_props["density"] = SingleAttrib('Density', 'val', [], [], [], missing_val=None)
OutputFile.output_props["MSD"] = SingleAttrib('MSD/Species', 'diffusionCoeff', [], [], ['-LMSD'], missing_val=None, skip_missing=True)
OutputFile.output_props["VACF"] = VACFOutputProperty()
OutputFile.output_props["NeventsSO"] = SingleAttrib('EventCounters/Entry[@Name="SOCells"]', # Outputfile tag name
                                                    'Count', # Outputfile tag attribute name
                                                    ["Rso"], # Required state variable
                                                    [], # Required output variables
                                                    [], # Required output plugins
                                                    div_by_N=True, # Divide the count by N
                                                    div_by_t=True, # Also divide by t
                                                    missing_val=0) # If counter is missing, return 0

if __name__ == "__main__":
    pass
