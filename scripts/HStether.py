#!/usr/bin/python3
import numpy, pydynamo, datastat
from pydynamo import ET

import math

################################################################
###      DEFINE THE "STATE" VARIABLES TO BE SWEPT & RANGE
################################################################
#This is the list of state variables and their ranges

densities = set(list(numpy.arange(0.1, 1.4, 0.1))+list(numpy.arange(0.9,1.05,0.01)))
densities = list(map(lambda x : datastat.roundSF(x, 3), list(densities)))
densities.sort()
Rso = list(map(lambda x : datastat.roundSF(x, 3), list(numpy.arange(0.01, 1.0, 0.01))))
statevars = [
    ("N", list(map(lambda x: 4*x**3, [5, 10, 15]))),
    ('ndensity', densities),
    ("Rso", Rso),
    ("InitState", ["FCC"]),
]

#You must also tell pydynamo about any new state variable (N,
#ndensity, InitState are built-in). First make a function that can
#figure out what the state variable is from a given config file.
######## \/   \/   \/  Below is already in pydynamo.py
# def Rso_config(XMLconfig):
#     tag = XMLconfig.tree.find('.//Global[@Type="SOCells"]')
#     if tag is None:
#         return float('inf')
#     else:
#         return pydynamo.conv_to_14sf(float(tag.attrib['Diameter']) / 2)

#Then tell pydynamo about this.
######## \/   \/   \/  Below is already in pydynamo.py
# pydynamo.ConfigFile.config_props["Rso"] = (True, Rso_config)
#The True value above tells pydynamo this state variable can be reliably
#rebuilt from a configuration file. Some values (like the initial
#lattice) are not able to be rebuilt.

#Rebuilding of state variables is what allows you to change the swept
#parameters of a simulation without having to rerun anything. If you
#add a state varible, then the value of the existing simulations can
#be determined automagically.

################################################################
###      DEFINE ANY OUTPUT TO BE COLLECTED
################################################################
# Here we want to collect a single value out from the output data
# file. Its attribute "Count" inside //EventCounters/Entry with the
# Name attribute value of "SOCells". We want this divided by the
# particle count AND the total run time.
######## \/   \/   \/  Below is already in pydynamo.py
# pydynamo.OutputFile.output_props["NeventsSO"] = pydynamo.SingleAttrib('EventCounters/Entry[@Name="SOCells"]', # Outputfile tag name
#                                                                      'Count', # Outputfile tag attribute name
#                                                                      ["Rso"], # Required state variable
#                                                                      [], # Required output variables
#                                                                      div_by_N=True, # Divide the count by N
#                                                                      div_by_t=True, # Also divide by t
#                                                                      missing_val=0) # If counter is missing, return 0
# pydynamo will handle weighted averaging this value from all existing
# simulations together, including propogating error estimates from
# jacknife averaging of different run sizes!

################################################################
###            DEFINE THE SETUP WORKER FUNCTION
################################################################
# This function is used (in parallel) to set up simulation configs for
# each state point.
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
        minR = max(0, (2**(2.5)*state['ndensity'])**(-1/3.0) - 0.5)
        #phiT= state['ndensity'] * (4/3) * math.pi * minR**3
        #minRho = max(0, (2**(1/6.0)-(6*state['ndensity']*(4/3)*minR**3)**(1/3))**3)
        if state['Rso'] <= minR:
            raise pydynamo.SkipThisPoint()

    # Again, for tethered systems in FCC lattices we do not simulate
    # much beyond a multiple of the minimum tether radius.
    if ("Rso" in state) and (state['ndensity'] >= 1.0) and (state['InitState'] == "FCC"):
        minR = max(0, (2**(2.5)*state['ndensity'])**(-1/3.0) - 0.5)
        if state['Rso'] >= 10*minR:
            raise pydynamo.SkipThisPoint()
        
    check_call(('dynamod -o '+config+' -m 0'+' -d ' + repr(state['ndensity'])+' -C'+str(Ncells)).split(), stdout=logfile, stderr=logfile)

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
mgr = pydynamo.SimManager("HSTetherWD", #Which subdirectory to work in
                          statevars, #State variables
                          ["p", "NeventsSO", "VACF"], # Output properties
                          restarts=3, #How many restarts (new initial configurations) should be done per state point
                          processes=None, #None is automatically use all processors
)

################################################################
###          REORGANISE ANY EXISTING SIMULATIONS
################################################################
# Just in case there's already data in the working dir, we smartly
# tidy it up depending on what state variables we're currently
# interested in.
#mgr.reorg_dirs()

################################################################
###          RUN SOME SIMULATIONS
################################################################
# Its safe to run this commmand as many times as you like, even if
# you're just processing data. It will just run more simulations to
# satisfy any missing values of the state variables and/or to bring
# the event count up to the required amount
mgr.run(setup_worker=setup_worker,
        particle_equil_events = 1000, # How many events per particle to equilibrate each sim for
        particle_run_events = 10000, # How many events per particle to run IN TOTAL
        particle_run_events_block_size=1000) # How big a block each run should be (for jacknife averaging).

# Note, changes to the block size only matters for simulations to be
# run, and will never cause run_events to be exceeded. Changes to
# equil_events only matters for state points without any simulations
# yet.

################################################################
###          GET THE DATA
################################################################
data = mgr.fetch_data() #This is a pandas dataframe with columns for
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
