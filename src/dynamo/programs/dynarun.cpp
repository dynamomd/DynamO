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
/*! \file dynarun.cpp 
 
  \brief Contains the main() function for dynarun
 
  Although this contains the main() function, most of the behaviour peculiar
  to dynarun is carried out by the Coordinator class.
*/

#include <dynamo/coordinator/coordinator.hpp>
#include <magnet/stream/formattedostream.hpp>
#include <magnet/arg_share.hpp>
#include <boost/program_options.hpp>
#include <iostream>
#include <signal.h>

/*! \brief Starting point for the dynarun program.
 
  This merely boots the Coordinator class. The Coordinator class is
  the true "main" function for the simulations.
 
  \param argc The number of command line arguments.
  \param argv A pointer to the array of command line arguments.
*/
int main(int argc, char *argv[])
{
  //Output the program licence
  std::cout << "dynarun  Copyright (C) 2011  Marcus N Campbell Bannerman\n"
	    << "This program comes with ABSOLUTELY NO WARRANTY.\n"
	    << "This is free software, and you are welcome to redistribute it\n"
	    << "under certain conditions. See the licence you obtained with\n"
	    << "the code\n";

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

      dynamo::Coordinator::get().parseOptions(args,argv);
      dynamo::Coordinator::get().initialise();
      dynamo::Coordinator::get().runSimulation();
      dynamo::Coordinator::get().outputData();
      dynamo::Coordinator::get().outputConfigs();
      
      return 0;
    }
  catch (std::exception& cep)
    {
      std::cout.flush();
      magnet::stream::FormattedOStream os(magnet::console::bold()
					  + magnet::console::red_fg() 
					  + "Main(): " + magnet::console::reset(), std::cerr);
      os << cep.what() << std::endl;
#ifndef DYNAMO_DEBUG
      os << "Try using the debugging executable for more information on the error." << std::endl;
#endif
      return 1;
    }
}
