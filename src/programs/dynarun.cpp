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
/*! \file mdrun.cpp 
 *
 * \brief Contains the main() function for dynarun
 *
 * Although this contains the main() function, most of the behaviour peculiar
 * to dynarun is carried out by the Coordinator class.
 */

#include "../src/coordinator/coordinator.hpp"
#include <buildinfo.hpp>
#include <magnet/arg_share.hpp>
#include <boost/program_options.hpp>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <signal.h>

/*! \brief Starting point for the dynarun program.
 *
 * This merely boots the Coordinator class. The Coordinator class is
 * the true "main" function for the simulations.
 *
 * \param argc The number of command line arguments.
 * \param argv A pointer to the array of command line arguments.
 */
int
main(int argc, char *argv[])
{
  /*! \brief The programs single instantiation of the simulation control class.
   */
  Coordinator coord;

  //Output the program licence
  std::cout << "dynarun  Copyright (C) 2011  Marcus N Campbell Bannerman\n"
	    << "This program comes with ABSOLUTELY NO WARRANTY.\n"
	    << "This is free software, and you are welcome to redistribute it\n"
	    << "under certain conditions. See the licence you obtained with\n"
	    << "the code\n"
	       "Git Checkout Hash " << GITHASH << "\n\n";

  //Reasonable precision for periodic output
  std::cout << std::setprecision(std::numeric_limits<float>::digits10);

  //Run the simulation
  try 
    {      
      magnet::ArgShare::getInstance().setArgs(argc, argv);
      
      //Strip off all arguments after our special -GLGTK marker
      
      int args = argc;
      for (int i(0); i < argc; ++i)
	if (!strcmp("-GLGTK", argv[i]))
	  {
	    args = i;
	    break;
	  }

      coord.parseOptions(args,argv);

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
