#!/usr/bin/python3
import numpy, pydynamo, datastat
from pydynamo import ET

import math

################################################################
###      DEFINE THE "STATE" VARIABLES TO BE SWEPT & RANGE
################################################################
#This is the list of state variables and their ranges

d=2

densities = set(list(numpy.arange(0.1, 0.9, 0.1))+list(numpy.arange(0.8,0.95,0.01)))
densities = list(map(lambda x : datastat.roundSF(x, 3), list(densities)))
densities.sort()
Rso = list(map(lambda x : datastat.roundSF(x, 3), list(numpy.arange(0.01, 2.0, 0.02))))

statevars = [
    ("N", list(map(lambda x: x**2, [12, 50, 100]))),
    ('ndensity', densities),
    ("Rso", Rso),
    ("InitState", ["SC", "hexagonal"]),
]

def setup_worker( config, #The name of the config file to generate.
                  state, #A dictionary of state variables to use
                  logfile, #File handle where to write progress/logging output
                  particle_equil_events, # How many events will be run per particle to equilibrate the config. Useful if in setup you also need to equilibrate an intermediate configuration.
):
    from subprocess import check_call
    
    Ncells_unrounded = state['N'] ** (1.0 / d)
    Ncells = int(round(Ncells_unrounded))
    if abs(Ncells - Ncells_unrounded) > 0.1:
        raise RuntimeError("Could not make "+str(state['N'])+" particles in an "+state['InitState']+" packing")

    check_call(('dynamod -m 0 -z 1 --rectangular-box --zero-vel 2  --i1 2 -T 1.0 -o '+config+' -d ' + repr(state['ndensity'])+' -x '+str(Ncells)+' -y '+str(Ncells)).split(), stdout=logfile, stderr=logfile)

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

    if state['InitState'] == "hexagonal":
        xml = pydynamo.ConfigFile(config)
        particles = xml.tree.findall(".//Pt")
        dx = (float(particles[1].find('P').attrib['x']) - float(particles[0].find('P').attrib['x'])) / 2
        for idx,pt in enumerate(particles):
            if ((idx // Ncells) % 2) == 1:
                p = pt.find('P')
                p.attrib['x']  = repr(float(p.attrib['x']) + dx)
        xml.save(config)
        

################################################################
###          CREATE A SIMULATION MANAGER
################################################################
mgr = pydynamo.SimManager("HDTetherWD", #Which subdirectory to work in
                          statevars, #State variables
                          ["p", "NeventsSO", "VACF"], # Output properties
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
mgr.run(setup_worker=setup_worker,
        particle_equil_events = 1000, # How many events per particle to equilibrate each sim for
        particle_run_events = 10000, # How many events per particle to run IN TOTAL
        particle_run_events_block_size=1000) # How big a block each run should be (for jacknife averaging).

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
