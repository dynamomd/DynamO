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
/*! \file mdrun.cpp 
 *
 * \brief Contains the main() function for dynarun
 *
 * Although this contains the main() function, most of the behaviour peculiar
 * to dynarun is carried out by the Coordinator class.
 */

#include <iostream>
#include <iomanip>
#include <signal.h>
#include <boost/program_options.hpp>
#include "../src/coordinator/coordinator.hpp"

/*! \brief The programs single instantiation of the simulation control class.
 */
Coordinator coord;


/*! \brief A function that merely wraps the Coordinator::signal_handler.

  \param sigtype The type of signal that has been recieved.
 */
void sig_handler_helper(int sigtype)
{
  coord.signal_handler(sigtype);
}

/*! \brief Starting point for the dynarun program.
 *
 * This merely registers some signal handlers and boots the
 * Coordinator class. The Coordinator class is the true "main"
 * function for the simulations.
 *
 * \param argc The number of command line arguments.
 * \param argv A pointer to the array of command line arguments.
 */
int
main(int argc, char *argv[])
{
  //Output the program licence
  std::cout << "dynarun  Copyright (C) 2009  Marcus N Campbell Bannerman\n"
	    << "This program comes with ABSOLUTELY NO WARRANTY.\n"
	    << "This is free software, and you are welcome to redistribute it\n"
	    << "under certain conditions. See the licence you obtained with\n"
	    << "the code\n\n";

  //Reasonable precision for periodic output
  std::cout << std::setprecision(std::numeric_limits<float>::digits10);

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
  catch (std::exception& cep)
    {
      fflush(stdout);
      std::cerr << cep.what();
      std::cerr << "\n" << IC_red 
		<< "MAIN:" << IC_reset << " Reached Main Error Loop"
#ifndef DYNAMO_DEBUG
		<< IC_red << "\nMAIN:" << IC_reset << "If this error is vauge, try using the debugging executable"
#endif
		<< "\n";

      return 1;
    }
}
