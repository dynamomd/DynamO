/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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
/*! \file coordinator.cpp 
 *
 * \brief Contains the code for the CCoordinator class.
 */

#include "coordinator.hpp"
#include "engine/include.hpp"
#include <signal.h>
#ifdef DYNAMO_DEBUG
# include <typeinfo>
#endif

void 
CCoordinator::signal_handler(int sigtype)
{
  switch (sigtype)
    {
    case SIGUSR1:
      //About to be stopped, fine by me, no mpi etc concerns yet
      return;
    case SIGUSR2:
      //Just try and shutdown before we're (kill -9)ed
      Engine->forceShutdown();
      return;
    case SIGINT:
      {
	//Clear the writes to screen
	std::cout << std::flush;
	std::cerr << "\n<S>hutdown, <E>xit, <D>ata or <P>eek at data output:";

	char c;
	//Clear the input buffer
	std::cin >> c;
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	switch (c)
	  {
	  case 's':
	  case 'S':
	    Engine->forceShutdown();
	    break;
	  case 'e':
	  case 'E':
	    exit(1);
	    break;
	  case 'p':
	  case 'P':
	    Engine->peekData();
	    break;
	  case 'd':
	  case 'D':
	    Engine->printStatus();
	    break;
	  }
	break;
      }
    }
}

boost::program_options::variables_map& 
CCoordinator::parseOptions(int argc, char *argv[])
{
  namespace po = boost::program_options;

  boost::program_options::options_description allopts(""),
    basicOpts, detailedEngineOpts, systemopts("System Options"), 
    engineopts("Engine Options")
    ;
  
  systemopts.add_options()
    ("help", "Produces this message")   
    ("n-threads,N", po::value<unsigned int>()->default_value(0), 
     "Number of threads to spawn for concurrent processing")
    ("out-config-file,o", po::value<std::string>(), 
     "Default config output file,(config.%ID.end.xml.bz2)")
    ("out-data-file", po::value<std::string>(), 
     "Default result output file (output.%ID.xml.bz2)")
    ("config-file", po::value<std::vector<std::string> >(), 
     "Specify a config file to load, or just list them on the command line")
    ;

  engineopts.add_options()
    ("engine-help", "Detailed options for the available engines")
    ("engine",boost::program_options::value<size_t>()->default_value(1),
     "Select Engine for simulation:\n"
     " Values:\n"
     "  1: \tSingle System Engine\n"
     "  2: \tNVT Replica Exchange\n"
     "  3: \tCompression Dynamics")
    ;

  basicOpts.add(systemopts).add(engineopts);

  CEngine::getCommonOptions(detailedEngineOpts);
  CEReplexer::getOptions(detailedEngineOpts);
  CECompressor::getOptions(detailedEngineOpts);
  
  allopts.add(basicOpts).add(detailedEngineOpts);

  boost::program_options::positional_options_description p;
  p.add("config-file", -1);
  
  boost::program_options::store(po::command_line_parser(argc, argv).
	    options(allopts).positional(p).run(), vm);
  boost::program_options::notify(vm);

  if (vm.count("help") || (argc==1)) 
    {
      std::cout << "Usage : mdrun <OPTION>...<config-file(s)>\n"
		<< "Initialises a configuration or loads a previous "
		<< "configuration, then calculates or loads the " 
		<< "trajectory and outputs data\n"
		<< basicOpts << "\n";
      exit(1);
    }

  if (vm.count("engine-help")) 
    {
      std::cout << "Engine Options:-"
		<< detailedEngineOpts << "\n";
      exit(1);
    }
  
  if (vm.count("config-file") == 0)
    D_throw() << "No configuration files to load specified";

  return vm;
}

void 
CCoordinator::initialise()
{
  switch (vm["engine"].as<size_t>())
    {
    case (1):
      Engine.set_ptr(new CESingle(vm));
      break;
    case (2):
      Engine.set_ptr(new CEReplexer(vm));
      break;
    case (3):
      Engine.set_ptr(new CECompressor(vm));
      break;
    default:
      D_throw() <<"Unknown Engine Selected"; 
    }
  
  Engine->initialisation();
}

void 
CCoordinator::runSimulation()
{
  //Only Run if there are collisions to run
  if (vm["ncoll"].as<unsigned long long>())
    Engine->runSimulation();
}

void 
CCoordinator::outputData()
{
  Engine->outputData();
}

void 
CCoordinator::outputConfigs()
{
  Engine->finaliseRun();

  //Only output if there are collisions to run
  if (vm["ncoll"].as<unsigned long long>())
    Engine->outputConfigs();
}
