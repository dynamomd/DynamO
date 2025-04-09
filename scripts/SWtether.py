#!/usr/bin/env python3
import math

import numpy

import pydynamo


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
        "HCP":4,
        "BCC":2,
        "SC":1,
    }
    Ncells_unrounded = (state['N'] / unitcellN[state['InitState']]) ** (1.0 / 3.0)
    Ncells = int(round(Ncells_unrounded))
    if abs(Ncells - Ncells_unrounded) > 0.1:
        raise RuntimeError("Could not make "+str(state['N'])+" particles in an "+state['InitState']+" packing")

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


    # Thermostat
    options = ''
    if 'kT' in state:
        options = options + ' -T '+repr(state['kT'])

    # Crystal lattice packing
    packmode = {
        'FCC':0,
        'BCC':1,
        'SC': 2,
        'HCP':3,
    }
    options = options + ' --i1 '+str(packmode[state['InitState']]) +' -C '+str(Ncells)
    
    # Square well or hard sphere?
    if state['Lambda'] != float('inf'):
        options = options + ' -m 1 + --f1 '+repr(state['Lambda'])
    else:
        options = options + ' -m 0'

    # denisty
    options = options + ' -d ' + repr(state['ndensity'])

    # Execution of dynamod
    check_call(('dynamod'+options+' -o '+config).split(), stdout=logfile, stderr=logfile)

    # Run of an equilibration step to blur the system state
    if ('Rso' in state) and (state['Rso'] != float('inf')) and (state["InitState"] == "Liquid"):
        print("\n", file=logfile)
        print("################################", file=logfile)
        print("#      Liquifaction Run        #", file=logfile)
        print("################################\n", file=logfile, flush=True)
        check_call(["dynarun", "--unwrapped", config, '-o', config, '-c', str(state['N'] * particle_equil_events), "--out-data-file", "data.liqequil.xml.bz2"], stdout=logfile, stderr=logfile)
    
    # Add the SO Cells global interaction (if needed)
    if ('Rso' in state) and (state['Rso'] != float('inf')):
        xml = pydynamo.ConfigFile(config)
        XMLGlobals = xml.tree.find(".//Globals")
        XMLSOCells = pydynamo.ET.SubElement(XMLGlobals, 'Global')
        XMLSOCells.attrib['Name'] = "SOCells" #Name can be anything
        XMLSOCells.attrib['Type'] = "SOCells" #This must be the right type of Global to load
        XMLSOCellsRange = pydynamo.ET.SubElement(XMLSOCells, 'Range')
        XMLSOCellsRange.attrib["Type"] = "All"
        XMLSOCells.attrib['Diameter'] = str(2 * state['Rso'])
        xml.save(config)


################################################################
###      DEFINE THE "STATE" VARIABLES TO BE SWEPT & RANGE
################################################################
#This is the list of state variables and their ranges

prefix="SWTether2"

statevars = [
#    [ #Sweep 
#        ("Lambda", [2]),
#        ("InitState", ["FCC"]),
#        ("PhiT", [float('inf')]),
#        ("N", list(map(lambda x: 4*x**3, [10]))), #15
#        ('ndensity', list(set(map(lambda x : pydynamo.roundSF(x, 3), list(numpy.arange(1.0, 1.41, 0.01)))))),
#        ("kT", list(set(map(lambda x : pydynamo.roundSF(x, 3), list(numpy.arange(1.0, 3.1, 0.1)))))),
#    ],
#    [ #Isochore
#        ("Lambda", [2]),
#        ("InitState", ["FCC"]),
#        ("PhiT", [float('inf')]),
#        ("N", list(map(lambda x: 4*x**3, [10]))), #15
#        ('ndensity', [1.3]),
#        ("kT", list(set(map(lambda x : pydynamo.roundSF(1/x, 3), list(numpy.arange(0.01, 1.01, 0.01)))))),
#    ],
    [ #Isochore
        ("Lambda", [1.5, 1.57]),
        ("InitState", ["HCP"]), #'FCC'
        ("PhiT", [float('inf')]),
        ("N", list(map(lambda x: 4*x**3, [10]))), #15
        ('ndensity', [1.3]),
        ("kT", list(set(map(lambda x : pydynamo.roundSF(1/x, 3), list(numpy.arange(0.01, 1.01, 0.01)))))),
    ],
    [ #Sweep 
        ("Lambda", [1.5, 1.57]),
        ("InitState", ["HCP"]), #'FCC'
        ("PhiT", [float('inf')]),
        ("N", list(map(lambda x: 4*x**3, [10]))), #15
        ('ndensity', list(set(map(lambda x : pydynamo.roundSF(x, 3), list(numpy.arange(1.0, 1.41, 0.01)))))),
        ("kT", list(set(map(lambda x : pydynamo.roundSF(x, 3), list(numpy.arange(0.5, 3.1, 0.1)))))),
    ],
#    [ #liquid runs
#        ("Lambda", [1.5, 1.57]),
#        ("InitState", ["HCP"]), #'FCC'
#        ("PhiT", [float('inf')]),
#        ("N", list(map(lambda x: 4*x**3, [10]))), #15
#        ('ndensity', list(set(map(lambda x : pydynamo.roundSF(x, 3), list(numpy.arange(0.01, 1.41, 0.01)))))),
#        ("kT", [5,4,3,2.5,1.5]),
#    ],
]
        
################################################################
###          CREATE A SIMULATION MANAGER
################################################################
mgr = pydynamo.SimManager("SWTether2", #Which subdirectory to work in
                          statevars, #State variables
                          ["p", "NeventsSO", 'cv', 'u',], # 'RadialDist' "VACF",  # Output properties
                          restarts=2, #How many restarts (new initial configurations) should be done per state point
                          processes=None, #None is automatically use all processors
)

################################################################
###          REORGANISE ANY EXISTING SIMULATIONS
################################################################
#mgr.reorg_dirs()

################################################################
###          RUN SOME SIMULATIONS
################################################################
#mgr.run(setup_worker=setup_worker,
#        particle_equil_events = 1000, # How many events per particle to equilibrate each sim for
#        particle_run_events = 10000, # How many events per particle to run IN TOTAL
#        particle_run_events_block_size=1000) # How big a block each run should be (for jacknife averaging).

################################################################
###          GET THE DATA
################################################################
#This creates a pandas dataframe with columns for the state variables
#AND any output values. It also generates pkl files, some for
#different properties.
data = mgr.fetch_data(1000)


