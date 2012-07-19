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

#include <dynamo/coordinator/engine/single.hpp>

namespace dynamo {
  ESingleSimulation::ESingleSimulation(const boost::program_options::variables_map& nVM, 
				       magnet::thread::ThreadPool& tp):
    Engine(nVM, "config.out.xml.bz2", "output.xml.bz2", tp),
    peekMode(false)
  {}

  void 
  ESingleSimulation::peekData()
  {
    peekMode = true;
    simulation.simShutdown();
  }

  void
  ESingleSimulation::runSimulation()
  {
    try {
      do
	{
	  if (peekMode)
	    {
	      simulation.endEventCount = vm["events"].as<size_t>();
	      simulation.outputData("peek.data.xml.bz2");
	      peekMode = false;
	    }
	
	  simulation.runSimulation();
	}
      while (peekMode);
    }
    catch (std::exception& cep)
      {
	try {
	  std::cerr << "\nEngine: Trying to output config to config.error.xml.bz2";
	  simulation.writeXMLfile("config.error.xml.bz2", !vm.count("unwrapped"));
	} catch (...)
	  {
	    std::cerr << "\nEngine: Could not output error config";
	  }
	throw;
      }
  }

  void
  ESingleSimulation::initialisation()
  {
    preSimInit();

    if (!(vm.count("config-file")) || 
	(vm["config-file"].as<std::vector<std::string> >().size() != 1))
      M_throw() << "You must only provide one input file in single mode";

    setupSim(simulation, vm["config-file"].as<std::vector<std::string> >()[0]);

    simulation.initialise();

    postSimInit(simulation);

    if (vm.count("ticker-period"))
      simulation.setTickerPeriod(vm["ticker-period"].as<double>());

  }

  void
  ESingleSimulation::outputData()
  {
    simulation.outputData(outputFormat.c_str());
  }

  void
  ESingleSimulation::outputConfigs()
  {
    simulation.writeXMLfile(configFormat.c_str(), !vm.count("unwrapped"));
  }

  void 
  ESingleSimulation::forceShutdown()
  {
    simulation.simShutdown();
  }
}
