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
/*! \file coordinator.cpp 
 *
 * \brief Contains the code for the Coordinator class.
 */

#include <dynamo/coordinator/coordinator.hpp>
#include <dynamo/coordinator/engine/include.hpp>
#include <cstdio>

namespace dynamo {
  Coordinator* Coordinator::_signal_handler = NULL;

  void 
  Coordinator::signal_handler(int sigtype)
  {
    //Restore the old method
    sigaction (SIGINT, &(_signal_handler->_old_SIGINT_handler), NULL);

    switch (sigtype)
      {
      case SIGUSR1:
	//About to be stopped, fine by me, no mpi etc concerns yet
	return;
      case SIGUSR2:
	//Just try and shutdown before we're (kill -9)ed
	_signal_handler->_engine->forceShutdown();
	return;
      case SIGINT:
	{
	  //Clear the writes to screen
	  std::cout.flush();
	  std::cerr << "\n<S>hutdown, <D>ata or <P>eek at data output:";

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
	      _signal_handler->_engine->forceShutdown();
	      break;
	      //	  case 'e':
	      //	  case 'E':
	      //	    if (_threads.getThreadCount())
	      //	      std::cerr << "Cannot <E>xit when threaded, causes program to hang. Try shutting down.";
	      //	    else
	      //	      exit(1);
	      //	    break;
	    case 'p':
	    case 'P':
	      _signal_handler->_engine->peekData();
	      break;
	    case 'd':
	    case 'D':
	      _signal_handler->_engine->printStatus();
	      break;
	    }
	  break;
	}
      }

    struct sigaction new_action;
    new_action.sa_handler = Coordinator::signal_handler;
    sigemptyset(&new_action.sa_mask);
    new_action.sa_flags = 0;
    sigaction (SIGINT, &new_action, NULL);
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
      ("n-threads,N", po::value<unsigned int>(),
       "Number of _threads to spawn for concurrent processing. (Only utilised by certain engine/sim configurations)")
      ("out-config-file,o", po::value<std::string>(),
       "Default config output file,(config.%ID.end.xml.bz2)")
      ("out-data-file", po::value<std::string>(),
       "Default result output file (output.%ID.xml.bz2)")
      ("config-file", po::value<std::vector<std::string> >(),
       "Specify a config file to load, or just list them on the command line")
      ;

    engineopts.add_options()
      ("engine",boost::program_options::value<size_t>()->default_value(1),
       "Select the Engine used to run the simulation:\n"
       " Values:\n"
       "  1: \tStandard Engine\n"
       "  2: \tNVT Replica Exchange Engine\n"
       "  3: \tCompression Engine")
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
	std::cout << "Usage : dynarun <OPTION>...<config-file(s)>\n"
		  << "Loads a configuration file, calculates the trajectory of the system using the specified engine and" 
		  << "outputs any collected data, including a the final configuration file.\n"
		  << basicOpts << "\n"
		  << detailedEngineOpts << "\n";
	exit(1);
      }

  
    if (vm.count("config-file") == 0)
      M_throw() << "No configuration files to load specified";

    return vm;
  }

  void 
  Coordinator::initialise()
  {
    //Register the signal handlers so we can respond to
    //attempts/warnings that the program will be killed
    {
      if (_signal_handler != NULL)
	M_throw() << "Can only have one instance of the Coordinator!";
    
      _signal_handler = this;
    
      //Build the handler response
      struct sigaction new_action;
      new_action.sa_handler = Coordinator::signal_handler;
      sigemptyset(&new_action.sa_mask);
      new_action.sa_flags = 0;
    
      //This is for Ctrl-c events
      sigaction (SIGINT, NULL, &_old_SIGINT_handler);
      if (_old_SIGINT_handler.sa_handler != SIG_IGN)
	sigaction (SIGINT, &new_action, NULL);
    
      struct sigaction old_action;
      //Sun Grid Engine sends this before a SIGSTOP if -notify is passed
      //to qsub
      sigaction (SIGUSR1, NULL, &old_action);
      if (old_action.sa_handler != SIG_IGN)
	sigaction (SIGUSR1, &new_action, NULL);
    
      //Sun Grid Engine sends this before a SIGKILL if -notify is passed
      //to qsub
      sigaction (SIGUSR2, NULL, &old_action);
      if (old_action.sa_handler != SIG_IGN)
	sigaction (SIGUSR2, &new_action, NULL);
    }

    if (vm.count("n-_threads"))
      _threads.setThreadCount(vm["n-_threads"].as<unsigned int>());

    switch (vm["engine"].as<size_t>())
      {
      case (1):
	_engine = shared_ptr<ESingleSimulation>(new ESingleSimulation(vm, _threads));
	break;
      case (2):
	_engine = shared_ptr<EReplicaExchangeSimulation>(new EReplicaExchangeSimulation(vm, _threads));
	break;
      case (3):
	_engine = shared_ptr<ECompressingSimulation>(new ECompressingSimulation(vm, _threads));
	break;
      default:
	M_throw() << vm["engine"].as<size_t>()
		  <<", Unknown Engine Number Selected"; 
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
}
