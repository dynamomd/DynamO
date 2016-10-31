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
/*! \file coordinator.cpp 
 *
 * \brief Contains the code for the Coordinator class.
 */

#include <dynamo/coordinator/coordinator.hpp>
#include <dynamo/coordinator/engine/include.hpp>
#include <cstdio>

//Need special treatment for the ways signals are handled on different platforms
#ifndef _WIN32
# include <signal.h>
#endif

namespace dynamo {
#ifdef _WIN32
  void
  Coordinator::setup_signal_handler() {
    if (!SetConsoleCtrlHandler((PHANDLER_ROUTINE) Coordinator::signal_handler, TRUE))
      M_throw() << "Failed to set signal control";
  }

  BOOL Coordinator::signal_handler(DWORD fdwCtrlType) {
    switch( fdwCtrlType ) 
      { 
	// Handle the CTRL-C signal. 
      case CTRL_C_EVENT:
	Coordinator::get()._engine->sigint();
	std::cerr << "\nCaught SIGINT, notifying running simulation...\n";
	return (TRUE);
      default:
	return FALSE;
      }
  }

#else
  void
  Coordinator::setup_signal_handler() {
    //Register the signal handlers so we can respond to
    //attempts/warnings that the program will be killed
      //Build the new handler response
      struct sigaction new_action;
      new_action.sa_handler = Coordinator::signal_handler;
      sigemptyset(&new_action.sa_mask);
      new_action.sa_flags = 0;
    
      //This is for Ctrl-c events
      struct sigaction old_action;
      sigaction (SIGINT, NULL, &old_action);
      if (old_action.sa_handler != SIG_IGN)
	sigaction (SIGINT, &new_action, NULL);

      sigaction (SIGTERM, NULL, &old_action);
      if (old_action.sa_handler != SIG_IGN)
	sigaction (SIGTERM, &new_action, NULL);
    }
    
  void 
  Coordinator::signal_handler(int sigtype)
  {
    switch (sigtype)
      {
      case SIGINT:
	{
	  //Disable this signal handler for any further SIGINT's, to
	  //let people kill the program with a double ctrl-c.
	  struct sigaction default_action;
	  default_action.sa_handler = SIG_DFL;
	  sigemptyset(&default_action.sa_mask);
	  default_action.sa_flags = SA_RESETHAND;
	  sigaction(SIGINT, &default_action, NULL);
	}
	Coordinator::get()._engine->sigint();
	std::cerr << "\nCaught SIGINT, notifying running simulation...\n";
	break;
      case SIGTERM:
	std::cerr << "\nCaught SIGTERM, notifying running simulation...\n";
	Coordinator::get()._engine->sigterm();
	break;
      }
  }
#endif

  boost::program_options::variables_map& 
  Coordinator::parseOptions(int argc, char *argv[])
  {
    namespace po = boost::program_options;

    boost::program_options::options_description allopts(""),
      basicOpts, detailedEngineOpts, systemopts("System Options"),
      engineopts("Engine Options")
      ;
  
#ifdef DYNAMO_bzip2_support
      std::string extension(".bz2");
#else
      std::string extension("");
#endif

      systemopts.add_options()
      ("help", "Produces this message")
      ("n-threads,N", po::value<unsigned int>(),
       "Number of threads to spawn for concurrent processing. (Only utilised by certain engine/sim configurations)")
      ("out-config-file,o", po::value<std::string>(),
       ("Default config output file,(config.%ID.end.xml"+extension+")").c_str())
      ("out-data-file", po::value<std::string>(),
       ("Default result output file (output.%ID.xml"+extension+")").c_str())
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
		  << "Loads a configuration file, calculates the dynamics of the system using the specified engine and" 
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
    setup_signal_handler();
    
    if (vm.count("n-threads"))
      _threads.setThreadCount(vm["n-threads"].as<unsigned int>());

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
  
    if (_enableVisualisation)
      _engine->enableVisualisation();
    _engine->initialisation();
  }

  void 
  Coordinator::runSimulation()
  {
    //Only run if there are collisions to run
    if (vm["events"].as<size_t>()) _engine->runSimulation();
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
    if (vm["events"].as<size_t>()) _engine->outputConfigs();
  }
}
