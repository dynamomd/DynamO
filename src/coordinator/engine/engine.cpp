/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "engine.hpp"
#include <limits>
#include "../../inputplugins/compression.hpp"
#include "../../dynamics/systems/tHalt.hpp"
#include "../../dynamics/systems/schedMaintainer.hpp"
#include "../../outputplugins/0partproperty/misc.hpp"
#include "../../outputplugins/general/reverseEvents.hpp"

void
CEngine::getCommonOptions(boost::program_options::options_description& opts)
{
  boost::program_options::options_description simopts("Common Engine Options");

  simopts.add_options()
    ("ncoll,c", boost::program_options::value<unsigned long long>()
     ->default_value(std::numeric_limits<unsigned long long>::max()),
     "No. of collisions in a trajectory")
    ("print-coll,p", boost::program_options::value<unsigned long long>()->default_value(100000), 
     "Default No. of collisions between periodic screen output")
    ("random-seed,s", boost::program_options::value<unsigned int>(),
     "Random seed for generator (To make the simulation reproduceable - Not for production use!)")
    ("ticker-period,t",boost::program_options::value<Iflt>(), 
     "Time between data collections. Defaults to the system MFT or 1 if no MFT available")
    ("scale-ticker",boost::program_options::value<Iflt>(), 
     "Useful when MFT data is available, can slow down or speed up the ticker in replex mode")
    ("equilibrate,E", "Turns off most output for a fast silent run")
    ("plugin-file,P", boost::program_options::value<std::string>(), "A list of output plugins to load")
    ("load-plugin,L", boost::program_options::value<std::vector<std::string> >(), "Additional individual plugins to load")
    ("halt-time,h", boost::program_options::value<Iflt>(),"Halt the system at this time")
    ("scheduler-maintainance,m", boost::program_options::value<Iflt>(),"Rebuild the scheduler"
     " periodically, for systems where we've not built the scheduler correctly")
    ;
  
  opts.add(simopts);
}


CEngine::CEngine(const boost::program_options::variables_map& nvm, 
		 std::string configFile, std::string outputFile,
		 CThreadPool& tp):
  vm(nvm),
  configFormat(configFile),
  outputFormat(outputFile),
  threads(tp)
{}

void 
CEngine::preSimInit()
{
  if (vm.count("out-config-file"))
    configFormat = vm["out-config-file"].as<std::string>();
  
  if (vm.count("out-data-file"))
    outputFormat = vm["out-data-file"].as<std::string>();
}


void 
CEngine::setupSim(CSimulation& Sim, const std::string filename)
{
  if (vm.count("random-seed"))
    Sim.setRandSeed(vm["random-seed"].as<unsigned int>());
  
  ////////////////////////Simulation Initialisation!!!!!!!!!!!!!
  //Now load the config
  Sim.loadXMLfile(filename.c_str());
  Sim.setTrajectoryLength(vm["ncoll"].as<unsigned long long>());
  
  if (vm["ncoll"].as<unsigned long long>() 
      > vm["print-coll"].as<unsigned long long>())
    Sim.setnPrint(vm["print-coll"].as<unsigned long long>());
  else
    Sim.setnPrint(vm["ncoll"].as<unsigned long long>());
    
  if (vm.count("halt-time"))
    Sim.addSystem(new CStHalt(&Sim, vm["halt-time"].as<Iflt>(), "SystemHaltEvent"));

  if (vm.count("scheduler-maintainance"))
    Sim.addSystem(new CSSchedMaintainer(&Sim, vm["scheduler-maintainance"].as<Iflt>(), "SchedulerRebuilder"));
  
  if (vm.count("plugin-file"))
    //Just add the plugins
    Sim.loadPlugins(vm["plugin-file"].as<std::string>());
  
  if (vm.count("load-plugin"))
    BOOST_FOREACH(const std::string& tmpString, vm["load-plugin"].as<std::vector<std::string> >())
      Sim.addOutputPlugin(tmpString);
  
  Sim.addOutputPlugin("ReverseEventsCheck");
  
  if (!vm.count("equilibrate"))
    //Just add the bare minimum outputplugin
    Sim.addOutputPlugin("Misc");
}
