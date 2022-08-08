#!/usr/bin/python3
import numpy, pydynamo, datastat
from pydynamo import ET

import math

def setup_worker( config, #The name of the config file to generate.
                  state, #A dictionary of state variables to use
                  logfile, #File handle where to write progress/logging output
                  particle_equil_events, # How many events will be run per particle to equilibrate the config. Useful if in setup you also need to equilibrate an intermediate configuration.
):
    from subprocess import check_call

    #Here we work out how many unit cells to make the system out of for various packings
    state = dict(state)
    if 'InitState' not in state:
        state['InitState'] = "FCC"
    
    unitcellN = {
        "FCC":4,
        "BCC":2,
        "SC":1,
    }
    Ncells_unrounded = (state['N'] / unitcellN[state['InitState']]) ** (1.0 / 3.0)
    Ncells = int(round(Ncells_unrounded))
    if abs(Ncells - Ncells_unrounded) > 0.1:
        raise RuntimeError("Could not make "+str(state['N'])+" particles in an "+state['InitState']+" packing")

    if state['ndensity'] <= 0 or ('kT' in state and state['kT'] <= 0):
        raise pydynamo.SkipThisPoint()
    
    # Here, for tethered systems, we do not simulate state points if
    # its going to be boring and "ideal". I only have worked out the
    # spacing expression for FCC, so all other crystals will just be
    # run regardless

    if ("Rso" in state) and (state['Rso'] != float('inf')) and (state['InitState'] == "FCC"):
        effrho = state['ndensity']*(state['Lambda']**3)
        minR = max(0, (2**(2.5)*effrho)**(-1/3.0) - 0.5)
        #phiT= state['ndensity'] * (4/3) * math.pi * minR**3
        #minRho = max(0, (2**(1/6.0)-(6*state['ndensity']*(4/3)*minR**3)**(1/3))**3)
        if state['Rso'] <= minR:
            raise pydynamo.SkipThisPoint()

    # This check is halting systems deep in the solid region, which should not be done!
    #
    #
    ### Again, for tethered systems in FCC lattices we do not simulate
    ### much beyond a multiple of the minimum tether radius.
    ##if ("Rso" in state) and (state['ndensity'] >= 1.0) and (state['InitState'] == "FCC"):
    ##    minR = max(0, (2**(2.5)*state['ndensity'])**(-1/3.0) - 0.5)
    ##    if state['Rso'] >= 10*minR:
    ##        raise pydynamo.SkipThisPoint()

    kTadd = ''
    if 'kT' in state:
        kTadd = '-T '+repr(state['kT'])
        
    check_call(('dynamod '+kTadd+' -o '+config+' -m 1'+' -d ' + repr(state['ndensity'])+' -C'+str(Ncells)+' --f1 '+repr(state['Lambda'])).split(), stdout=logfile, stderr=logfile)

    if ('Rso' in state) and (state['Rso'] != float('inf')) and (state["InitState"] == "Liquid"):
        print("\n", file=logfile)
        print("################################", file=logfile)
        print("#      Liquifaction Run        #", file=logfile)
        print("################################\n", file=logfile, flush=True)
        check_call(["dynarun", "--unwrapped", config, '-o', config, '-c', str(state['N'] * particle_equil_events), "--out-data-file", "data.liqequil.xml.bz2"], stdout=logfile, stderr=logfile)
    
    #Add the SO Cells global interaction (if needed)
    if ('Rso' in state) and (state['Rso'] != float('inf')):
        xml = pydynamo.ConfigFile(config)
        XMLGlobals = xml.tree.find(".//Globals")
        XMLSOCells = ET.SubElement(XMLGlobals, 'Global')
        XMLSOCells.attrib['Name'] = "SOCells" #Name can be anything
        XMLSOCells.attrib['Type'] = "SOCells" #This must be the right type of Global to load
        XMLSOCellsRange = ET.SubElement(XMLSOCells, 'Range')
        XMLSOCellsRange.attrib["Type"] = "All"
        XMLSOCells.attrib['Diameter'] = str(2 * state['Rso'])
        xml.save(config)


################################################################
###      DEFINE THE "STATE" VARIABLES TO BE SWEPT & RANGE
################################################################
#This is the list of state variables and their ranges

prefix="SW_states"

import pandas
df = pandas.read_csv("states.csv")
df = df.drop(columns=['Unnamed: 0']).rename(columns={'lam':'Lambda', 'T':'kT', 'density':'ndensity'})
df['InitState'] = "FCC"
df['N'] = 4 * 10**3
df = df[df['ndensity'] > 0]
states= df.to_dict('records')

################################################################
###          CREATE A SIMULATION MANAGER
################################################################
mgr = pydynamo.SimManager(workdir=prefix, #Which subdirectory to work in
                          statevars=None, #State variables
                          outputs=["p", 'cv', 'u'], # Output properties
                          restarts=2, #How many restarts (new initial configurations) should be done per state point
                          processes=None, #None is automatically use all processors
                          states=states
)

################################################################
###          REORGANISE ANY EXISTING SIMULATIONS
################################################################
#mgr.reorg_dirs()

################################################################
###          RUN SOME SIMULATIONS
################################################################
mgr.run(setup_worker=setup_worker,
        particle_equil_events = 1000, # How many events per particle to equilibrate each sim for
        particle_run_events = 3000, # How many events per particle to run IN TOTAL
        particle_run_events_block_size=1000) # How big a block each run should be (for jacknife averaging).

################################################################
###          GET THE DATA
################################################################
#This creates a pandas dataframe with columns for the state variables
#AND any output values. It also generates pkl files, some for
#different properties.
data = mgr.fetch_data(1000)


