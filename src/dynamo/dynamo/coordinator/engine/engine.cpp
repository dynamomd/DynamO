/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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
#include <dynamo/systems/tHalt.hpp>
#include <dynamo/systems/visualizer.hpp>
#include <limits>


namespace dynamo {
  void
  Engine::getCommonOptions(boost::program_options::options_description& opts)
  {
    boost::program_options::options_description simopts("Common Engine Options");

    simopts.add_options()
      ("events,c", boost::program_options::value<size_t>()
       ->default_value(std::numeric_limits<size_t>::max(), "no-limit"),
       "No. of events to run the simulation for.")
      ("print-events,p", boost::program_options::value<size_t>()->default_value(100000), 
       "No. of events between periodic screen output.")
      ("random-seed,s", boost::program_options::value<unsigned int>(),
       "Random seed for generator (To make the simulation reproduceable - Only for debugging!)")
      ("ticker-period,t",boost::program_options::value<double>(), 
       "Time between data collections. Defaults to the system MFT or 1 if no MFT available")
      ("equilibrate,E", "Turns off most output for a fast silent run")
      ("load-plugin,L", boost::program_options::value<std::vector<std::string> >(), 
       "Additional individual plugins to load")
      ("sim-end-time,f", boost::program_options::value<double>()->default_value(std::numeric_limits<double>::max(), "no limit"), 
       "Simulation end time (Note, In replica exchange, each systems end time is scaled by"
       "(T_cold/T_i)^{1/2}, see replex-interval)")
      ("unwrapped", "Don't apply the boundary conditions of the system when writing out the particle positions.")
      ("snapshot", boost::program_options::value<double>(),
       "Sets the system time inbetween saving snapshots of the system.")
      ("snapshot-events", boost::program_options::value<size_t>(),
       "Sets the event count inbetween saving snapshots of the system.")
      ;
  
    opts.add(simopts);
  }


  Engine::Engine(const boost::program_options::variables_map& nvm, 
		 std::string configFile, std::string outputFile,
		 magnet::thread::ThreadPool& tp):
    vm(nvm),
    configFormat(configFile),
    outputFormat(outputFile),
    _SIGINT(false),
    _SIGTERM(false),
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
    Sim.ranGenerator.seed(std::random_device()());
    if (vm.count("random-seed"))
      Sim.ranGenerator.seed(vm["random-seed"].as<unsigned int>());
  
    ////////////////////////Simulation Initialisation!!!!!!!!!!!!!
    //Now load the config
    Sim.loadXMLfile(filename.c_str());
    
    Sim.endEventCount = vm["events"].as<size_t>();
  
    if (vm["events"].as<size_t>() 
	> vm["print-events"].as<size_t>())
      Sim.eventPrintInterval = vm["print-events"].as<size_t>();
    else
      Sim.eventPrintInterval = vm["events"].as<size_t>();
    
    if (vm.count("sim-end-time") && (dynamic_cast<const EReplicaExchangeSimulation*>(this) == NULL))
      Sim.systems.push_back(shared_ptr<System>(new SystHalt(&Sim, vm["sim-end-time"].as<double>(), "SystemStopEvent")));

#ifdef DYNAMO_visualizer
    Sim.systems.push_back(shared_ptr<System>(new SVisualizer(&Sim, filename, Sim.lastRunMFT)));
#endif

    if (vm.count("load-plugin"))
      {
	for (const std::string& tmpString : vm["load-plugin"].as<std::vector<std::string> >())
	  Sim.addOutputPlugin(tmpString);
      }
  
    if (!vm.count("equilibrate"))
      //Just add the bare minimum outputplugin
      Sim.addOutputPlugin("Misc");
  }
}
