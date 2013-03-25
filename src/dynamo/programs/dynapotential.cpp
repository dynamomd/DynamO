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

#include <dynamo/interactions/potentials/lennard_jones.hpp>
#include <boost/program_options.hpp>

/*! \brief Starting point for the dynapotential program.
 
  \param argc The number of command line arguments.
  \param argv A pointer to the array of command line arguments.
*/
int main(int argc, char *argv[])
{
  //Output the program licence
  std::cout << "dynapotential  Copyright (C) 2011  Marcus N Campbell Bannerman\n"
	    << "This program comes with ABSOLUTELY NO WARRANTY.\n"
	    << "This is free software, and you are welcome to redistribute it\n"
	    << "under certain conditions. See the licence you obtained with\n"
	    << "the code\n";
  try 
    {      
      namespace po = boost::program_options;
      
      boost::program_options::variables_map vm;
      boost::program_options::options_description options("Program Options");
      
      options.add_options()
	("help", "Produces this message")   
	("cutoff", po::value<double>()->default_value(3), "The cutoff radius of the potential")
	("attractive-steps", po::value<double>()->default_value(1), "The number of steps in the attractive part of the potential")
	("steps", po::value<size_t>()->default_value(10), "The number of steps to output data for")
	("volume", "Use the volume averaged energy algorithm")
	("left", "Use the left energy algorithm")
	("mid", "Use the midpoint energy algorithm")
	("right", "Use the right energy algorithm")
	("B2", "Use the B2 algorithm")
	("kT", po::value<double>()->default_value(1), "Set the temperature for the B2 algorithm")
	;

      boost::program_options::store(po::command_line_parser(argc, argv).
				    options(options).run(), vm);
      boost::program_options::notify(vm);
    
      if (vm.count("help")) 
	{
	  M_throw() << "Usage : dynahist_rw <OPTION>...<data-file(s)>\n"
		    << "Determines the weighting functions for the histograms\n"
		    << options << "\n";
	}
      
      using namespace dynamo;

      std::cout.precision(15);
 
      PotentialLennardJones::UMode U_mode;
      PotentialLennardJones::RMode R_mode = PotentialLennardJones::DELTAU;

      if (vm.count("mid")) U_mode = PotentialLennardJones::MIDPOINT;
      else if (vm.count("left")) U_mode = PotentialLennardJones::LEFT;
      else if (vm.count("right")) U_mode = PotentialLennardJones::RIGHT;
      else if (vm.count("volume")) U_mode = PotentialLennardJones::VOLUME;
      else if (vm.count("B2")) U_mode = PotentialLennardJones::SECONDVIRIAL;
      else 
	M_throw() << "Could not find which energy algorithm to use";

      PotentialLennardJones deltaU(1.0, 1.0, vm["cutoff"].as<double>(), U_mode, R_mode, vm["attractive-steps"].as<double>(), vm["kT"].as<double>());
      
      for (size_t i(0); i < vm["steps"].as<size_t>(); ++i)
	std::cout << deltaU[i].first << " " << deltaU[i].second << "\n";
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
