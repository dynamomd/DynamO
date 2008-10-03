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
#include <iostream>
#include <signal.h>
#include <boost/program_options.hpp>
#include "../src/coordinator/coordinator.hpp"

CCoordinator coord;

/*! \file mdrun.cpp 
 *
 * Contains the main() function for dynarun.cpp

 */
void sig_handler_helper(int i)
{
  coord.signal_handler(i);
}

/*! Just the starting point of the program
 *
 * This is the starting point for the dynarun program.  program. This
 * merely registers some signal handlers and boots the CCoordinator
 * class. The CCoordinator class is the true "main" function.
 *
 */

int
main(int argc, char *argv[])
{
  //Output the program licence
  std::cout << "dynarun  Copyright (C) 2008  Marcus N Campbell Bannerman\n"
	    << "This program comes with ABSOLUTELY NO WARRANTY.\n"
	    << "This is free software, and you are welcome to redistribute it\n"
	    << "under certain conditions. See the licence you obtained with\n"
	    << "the code\n\n";

  //Register the signal handlers so we can respond to
  //attempts/warnings that the program will be killed
  {
    //Build the handler response
    struct sigaction new_action, old_action;
    new_action.sa_handler = sig_handler_helper;
    sigemptyset (&new_action.sa_mask);
    new_action.sa_flags = 0;
    
    //This is for Ctrl-c events
    sigaction (SIGINT, NULL, &old_action);
    if (old_action.sa_handler != SIG_IGN)
      sigaction (SIGINT, &new_action, NULL);
    
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

  //Run the simulation
  try 
    {      
      coord.parseOptions(argc,argv);

      coord.initialise();
      
      coord.runSimulation();

      coord.outputData();

      coord.outputConfigs();

      std::cout << "\n";

      return 0;
    }
  catch (DYNAMO::Exception& cep)
    {
      fflush(stdout);
      std::cerr << cep.what();
      std::cerr << "\n" << IC_red 
		<< "MAIN:" << IC_reset << " Reached Main Error Loop"
		<< "\n";

      return 1;
    }
}
