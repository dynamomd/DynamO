#!/usr/bin/python3
import numpy, pydynamo, datastat
from pydynamo import ET

import math

################################################################
###      DEFINE THE "STATE" VARIABLES TO BE SWEPT & RANGE
################################################################
#This is the list of state variables and their ranges

#Old SWTether runs
densities = list(set(map(lambda x : datastat.roundSF(x, 3), list(numpy.arange(0.05, 1.4, 0.1)))))
phi_T =     [float('inf')] + list(set(map(lambda x : datastat.roundSF(x, 3), [0.001, 0.005, 0.01, 0.02, 0.03, 0.04]+list(numpy.arange(0.05, 3.0, 0.05)))))


#New SWTether2 runs, isotherm
phi_T = [float('inf')]
densities = list(set(map(lambda x : datastat.roundSF(x, 3), list(numpy.arange(0.01, 1.3, 0.01)))))
kT=[1.0]

#New SWTether2 runs, isochor
phi_T = [float('inf')]
densities = [1.29]
kT=list(set(map(lambda x : datastat.roundSF(1/x, 3), list(numpy.arange(0.01, 1.00, 0.01)))))

statevars = [
    ("N", list(map(lambda x: 4*x**3, [10]))), #15
    ('ndensity', densities),
    #("Rso", Rso),
    ("PhiT", phi_T),
    ("Lambda", [2.0]),
    ("kT", kT),
    ("InitState", ["FCC"]),
]


def setup_worker( config, #The name of the config file to generate.
                  state, #A dictionary of state variables to use
                  logfile, #File handle where to write progress/logging output
                  particle_equil_events, # How many events will be run per particle to equilibrate the config. Useful if in setup you also need to equilibrate an intermediate configuration.
):
    from subprocess import check_call

    #Here we work out how many unit cells to make the system out of for various packings
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
###          CREATE A SIMULATION MANAGER
################################################################
mgr = pydynamo.SimManager("SWTether2", #Which subdirectory to work in
                          statevars, #State variables
                          ["p", "NeventsSO", "VACF", 'cv', 'u'], # Output properties
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
data = mgr.fetch_data(1000) #This is a pandas dataframe with columns for
                        #the state variables AND any output values

# You can just write it out as a spreadsheet
#data.to_csv("output.csv")
# Or pickle it for later processing
import pickle
pickle.dump(data, open("tether.pkl", 'wb'))

#import matplotlib.pyplot as plt
#
#plt.errorbar(density_slice["Rso"], list(map(lambda x: x.nominal_value, density_slice["NeventsSO"])), yerr=list(map(lambda x: x.std_dev, density_slice["NeventsSO"])), fmt='x')
#plt.plot(density_slice["Rso"], 3 / density_slice["Rso"]/math.sqrt(2*math.pi))
#plt.show()

#density_slice.plot(x="Rso", y="NeventsSO")

# pydynamo also has some nice functions that propogate uncertainty
# through operations like numerical integration
#from datastat import simpson_integrate
#integral_of_NeventsSO = simpson_integrate(density_slice['Rso'].values[:-1], density_slice['NeventsSO'].values[:-1])
#print(integral_of_NeventsSO)
# The [:-1] is needed to integrate the reversed function
