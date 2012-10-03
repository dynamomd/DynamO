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

#include <signal.h>
#include <dynamo/coordinator/coordinator.hpp>
#include <dynamo/coordinator/engine/single.hpp>

namespace dynamo {
  ESingleSimulation::ESingleSimulation(const boost::program_options::variables_map& nVM, 
				       magnet::thread::ThreadPool& tp):
    Engine(nVM, "config.out.xml.bz2", "output.xml.bz2", tp)
  {}

  void
  ESingleSimulation::runSimulation()
  {
    try {
      while (true)
	{
	  if (!simulation.runSimulationStep()) break;
	  if (_SIGINT)
	    {
	      //Clear the writes to screen
	      std::cout.flush();
	      std::cerr << "\n<S>hutdown or <P>eek at data output:";
	      
	      char c;
	      //Clear the input buffer
	      std::cin.clear();
	      setvbuf(stdin, NULL, _IONBF, 0);
	      c=getchar();
	      setvbuf(stdin, NULL, _IOLBF, 0);
	      switch (c)
		{
		case 's':
		case 'S':
		  simulation.simShutdown();
		  break;
		case 'p':
		case 'P':
		  simulation.outputData("peek.data.xml.bz2");
		  break;
		}	      
	      {
		struct sigaction new_action;
		new_action.sa_handler = dynamo::Coordinator::signal_handler;
		sigemptyset(&new_action.sa_mask);
		new_action.sa_flags = 0;
		sigaction(SIGINT, &new_action, NULL);
	      }
	    }
	}
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
}
