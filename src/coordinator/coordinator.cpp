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
/*! \file coordinator.cpp 
 *
 * \brief Contains the code for the Coordinator class.
 */

#include "coordinator.hpp"
#include "engine/include.hpp"
#include <signal.h>
#ifdef DYNAMO_DEBUG
# include <typeinfo>
#endif

void 
Coordinator::signal_handler(int sigtype)
{
  switch (sigtype)
    {
    case SIGUSR1:
      //About to be stopped, fine by me, no mpi etc concerns yet
      return;
    case SIGUSR2:
      //Just try and shutdown before we're (kill -9)ed
      _engine->forceShutdown();
      return;
    case SIGINT:
      {
	//Clear the writes to screen
	std::cout.flush();
	std::cerr << "\n<S>hutdown, <E>xit, <D>ata or <P>eek at data output:";

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
	    _engine->forceShutdown();
	    break;
	  case 'e':
	  case 'E':
	    if (_threads.getThreadCount())
	      std::cerr << "Cannot <E>xit when threaded, causes program to hang. Try shutting down.";
	    else
	      exit(1);
	    break;
	  case 'p':
	  case 'P':
	    _engine->peekData();
	    break;
	  case 'd':
	  case 'D':
	    _engine->printStatus();
	    break;
	  }
	break;
      }
    }
}

boost::program_options::variables_map& 
Coordinator::parseOptions(int argc, char *argv[])
{
  namespace po = boost::program_options;

  boost::program_options::options_description allopts(""),
    basicOpts, detailedEngineOpts, systemopts("System Options"), 
    engineopts("Engine Options")
    ;
  
  systemopts.add_options()
    ("help", "Produces this message")   
    ("n-_threads,N", po::value<unsigned int>(), 
     "Number of _threads to spawn for concurrent processing. (Only utilised by some engine/sim configurations)")
    ("out-config-file,o", po::value<std::string>(), 
     "Default config output file,(config.%ID.end.xml.bz2)")
    ("out-data-file", po::value<std::string>(), 
     "Default result output file (output.%ID.xml.bz2)")
    ("config-file", po::value<std::vector<std::string> >(), 
     "Specify a config file to load, or just list them on the command line")
    ("uncompressed", "Output the XML config file without bzip compression, you have to specify out-data-file and out-config-file if you use this option")
    ;

  engineopts.add_options()
    ("engine-help", "Detailed options for the available engines")
    ("engine",boost::program_options::value<size_t>()->default_value(1),
     "Select Engine for simulation:\n"
     " Values:\n"
     "  1: \tSingle System Engine\n"
     "  2: \tNVT Replica Exchange\n"
     "  3: \tCompression dynamics")
    ;

  basicOpts.add(systemopts).add(engineopts);

  Engine::getCommonOptions(detailedEngineOpts);
  EReplicaExchangeSimulation::getOptions(detailedEngineOpts);
  ECompressingSimulation::getOptions(detailedEngineOpts);
  
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

  if (vm.count("uncompressed") && !vm.count("out-config-file"))
    D_throw() << "When using uncompressed output you must specify the output config file name";

  if (vm.count("out-config-file"))
    {
      std::string fileName(vm["out-config-file"].as<std::string>());
      if (vm.count("uncompressed") 
	  && (std::string(fileName.end()-4, fileName.end()) == ".bz2"))
	D_throw() << "You should not use a .bz2 extension for uncompressed"
	  " config files";
    }
  
  
  if (vm.count("uncompressed") && !vm.count("out-data-file"))
    D_throw() << "When using uncompressed output you must specify the output data file name";
  
  if (vm.count("out-data-file"))
    {
      std::string fileName(vm["out-data-file"].as<std::string>());
      if (vm.count("uncompressed") 
	  && (std::string(fileName.end()-4, fileName.end()) == ".bz2"))
	D_throw() << "You should not use a .bz2 extension for uncompressed"
	  " output files";
    }
  
  return vm;
}

void 
Coordinator::initialise()
{
  if (vm.count("n-_threads"))
    _threads.setThreadCount(vm["n-_threads"].as<unsigned int>());

  switch (vm["engine"].as<size_t>())
    {
    case (1):
      _engine.set_ptr(new ESingleSimulation(vm, _threads));
      break;
    case (2):
      _engine.set_ptr(new EReplicaExchangeSimulation(vm, _threads));
      break;
    case (3):
      _engine.set_ptr(new ECompressingSimulation(vm, _threads));
      break;
    default:
      D_throw() <<"Unknown Engine Selected"; 
    }
  
  _engine->initialisation();
}

void 
Coordinator::runSimulation()
{
  //Only Run if there are collisions to run
  if (vm["ncoll"].as<unsigned long long>())
    _engine->runSimulation();

}

void 
Coordinator::outputData()
{
  _engine->outputData();
}

void 
Coordinator::outputConfigs()
{
  _engine->finaliseRun();

  //Only output if there are collisions to run
  if (vm["ncoll"].as<unsigned long long>())
    _engine->outputConfigs();
}
