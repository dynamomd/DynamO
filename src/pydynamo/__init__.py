#!/usr/bin/env python3
"""A python module to simplify supervision and running of multiple DynamO simulations.

You define state variables which you'd like to vary during a parameter
sweep as well as output properties you'd like to collect. This module
then manages running the simulations, sorting them if the state
variables have changed/increased/etc, continuing the simulations if
they were aborted or the parameters were changed, as well as
calculating error estimates for the properties obtained at each state
point.
"""
import glob
import os
import subprocess
import time
from multiprocessing import Pool, cpu_count

import alive_progress

from pydynamo.config_files import ConfigFile, validate_configfile
from pydynamo.file_types import ET
from pydynamo.math import conv_to_14sf, print_to_14sf, roundSF
from pydynamo.output_properties import OutputFile, validate_outputfile
from pydynamo.weighted_types import KeyedArray, WeightedType


class SkipThisPoint(BaseException):
    pass

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
            
    def imap_unordered(self, func, iterable, chunksize=2):
        """A version of imap_unordered that works with multiprocessing.Pool"""
        # Chunksize halves overhead for tiny tasks, but also doesn't limit
        # parallelism when doing a few slow tasks, or on a system with many
        # cores
        if self.processes == 1:
            # If we only have one process, just use the builtin map
            for item in iterable:
                yield func(item)
            return
        else:
            pool = Pool(processes=self.processes)
            for result in pool.imap_unordered(func, iterable, chunksize):
                yield result


    def reorg_dirs(self):
        entries = os.listdir(self.workdir)
        print("Reorganising existing data directories...")
        n = len(entries)
        
        with alive_progress.alive_bar(n) as progress:
            #This is a parallel loop, returning items as they finish in arbitrary order
            for actions in self.imap_unordered(reorg_dir_worker, [(d, self) for d in entries]):
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

        import collections

        #We store the extracted data in a dict of dicts. The first
        #dict is for the state, the second for the property.
        state_data = collections.defaultdict(dict)        

        #So we run the per data dir operation, then reduce everything
        state_data = {}
        with alive_progress.alive_bar(n) as progress:
            #This is a parallel loop, returning items as they finish in arbitrary order
            for result in self.imap_unordered(perdir, [(d, particle_equil_events, self) for d in output_dirs]):
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
        
        import pickle
        pickle.dump(state_data, open(self.workdir+".raw_data.pkl", 'wb'))

        #We now prep the data for processing, we convert our
        #WeightedType's to ufloats as pandas supports that natively.
        for statevars, data in state_data.items():
            for prop in data:
                if isinstance(data[prop], WeightedType) and not isinstance(data[prop], KeyedArray):
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

        ##Now we write out the data

        return df
    