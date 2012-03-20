/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

    This program is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    version 3 as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
 
#include <dynamo/coordinator/engine/engine.hpp>
#include <dynamo/coordinator/engine/replexer.hpp>
#include <dynamo/inputplugins/compression.hpp>
#include <dynamo/dynamics/systems/tHalt.hpp>
#include <dynamo/dynamics/systems/schedMaintainer.hpp>
#include <dynamo/outputplugins/0partproperty/misc.hpp>
#include <dynamo/dynamics/systems/visualizer.hpp>
#include <dynamo/dynamics/systems/snapshot.hpp>
#include <limits>


namespace dynamo {
  void
  Engine::getCommonOptions(boost::program_options::options_description& opts)
  {
    boost::program_options::options_description simopts("Common Engine Options");

    simopts.add_options()
      ("events,c", boost::program_options::value<unsigned long long>()
       ->default_value(std::numeric_limits<unsigned long long>::max(), "no-limit"),
       "No. of events to run the simulation for.")
      ("print-events,p", boost::program_options::value<unsigned long long>()->default_value(100000), 
       "No. of events between periodic screen output.")
      ("random-seed,s", boost::program_options::value<unsigned int>(),
       "Random seed for generator (To make the simulation reproduceable - Only for debugging!)")
      ("ticker-period,t",boost::program_options::value<double>(), 
       "Time between data collections. Defaults to the system MFT or 1 if no MFT available")
#ifdef DYNAMO_visualizer    
      ("visualizer,V", 
       "Enables the visualizer and sets the initial update frequency")
#endif
      ("equilibrate,E", "Turns off most output for a fast silent run")
      ("load-plugin,L", boost::program_options::value<std::vector<std::string> >(), 
       "Additional individual plugins to load")
      ("sim-end-time,f", boost::program_options::value<double>()->default_value(std::numeric_limits<double>::max(), "no limit"), 
       "Simulation end time (Note, In replica exchange, each systems end time is scaled by"
       "(T_cold/T_i)^{1/2}, see replex-interval)")
      ("scheduler-maintainance,m", boost::program_options::value<double>(),
       "Rebuild the scheduler periodically, for systems where we've not built "
       "the scheduler correctly")
      ("unwrapped", "Don't apply the boundary conditions of the system when writing out the particle positions.")
      ("snapshot", boost::program_options::value<double>(),
       "Sets the system time inbetween saving snapshots of the system.")
      ;
  
    opts.add(simopts);
  }


  Engine::Engine(const boost::program_options::variables_map& nvm, 
		 std::string configFile, std::string outputFile,
		 magnet::thread::ThreadPool& tp):
    vm(nvm),
    configFormat(configFile),
    outputFormat(outputFile),
    threads(tp)
  {}

  void 
  Engine::preSimInit()
  {
    if (vm.count("out-config-file"))
      configFormat = vm["out-config-file"].as<std::string>();
  
    if (vm.count("out-data-file"))
      outputFormat = vm["out-data-file"].as<std::string>();
  }


  class EReplicaExchangeSimulation;

  void 
  Engine::setupSim(Simulation& Sim, const std::string filename)
  {
    if (vm.count("random-seed"))
      Sim.setRandSeed(vm["random-seed"].as<unsigned int>());
  
    ////////////////////////Simulation Initialisation!!!!!!!!!!!!!
    //Now load the config
    Sim.loadXMLfile(filename.c_str());
    Sim.configLoaded();
    Sim.setTrajectoryLength(vm["events"].as<unsigned long long>());
  
    if (vm["events"].as<unsigned long long>() 
	> vm["print-events"].as<unsigned long long>())
      Sim.setnPrint(vm["print-events"].as<unsigned long long>());
    else
      Sim.setnPrint(vm["events"].as<unsigned long long>());
    
    if (vm.count("sim-end-time") && (dynamic_cast<const EReplicaExchangeSimulation*>(this) == NULL))
      Sim.addSystem(shared_ptr<System>(new SystHalt(&Sim, vm["sim-end-time"].as<double>(), "SystemStopEvent")));

    if (vm.count("scheduler-maintainance"))
      Sim.addSystem(shared_ptr<System>(new SysSchedMaintainer(&Sim, vm["scheduler-maintainance"].as<double>(), "SchedulerRebuilder")));

#ifdef DYNAMO_visualizer
    if (vm.count("visualizer"))
      Sim.addSystem(shared_ptr<System>(new SVisualizer(&Sim, filename, Sim.lastRunMFT)));
#endif

    if (vm.count("snapshot"))
      Sim.addSystem(shared_ptr<System>(new SSnapshot(&Sim, vm["snapshot"].as<double>(), "SnapshotEvent")));

    if (vm.count("load-plugin"))
      {
	BOOST_FOREACH(const std::string& tmpString, 
		      vm["load-plugin"].as<std::vector<std::string> >())
	  Sim.addOutputPlugin(tmpString);
      }
  
    if (!vm.count("equilibrate"))
      //Just add the bare minimum outputplugin
      Sim.addOutputPlugin("Misc");
  }
}
