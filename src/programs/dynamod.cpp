/*  DYNAMO:- Event driven molecular dynamics simulator 
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
#include <iostream>
#include <cstdio>
#include <signal.h>
#include <boost/program_options.hpp>
using namespace std;
using namespace boost;
namespace po = boost::program_options;

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include "../src/simulation/simulation.hpp"
#include "../src/dynamics/BC/include.hpp"
#include "../src/dynamics/dynamics.hpp"
#include "../src/dynamics/systems/ghost.hpp"
#include <magnet/exception.hpp>
#include "../src/schedulers/include.hpp"
#include "../src/inputplugins/include.hpp"
#include <buildinfo.hpp>
Simulation sim;

int
main(int argc, char *argv[])
{
  std::cout << "dynamod  Copyright (C) 2011  Marcus N Campbell Bannerman\n"
	    << "This program comes with ABSOLUTELY NO WARRANTY.\n"
	    << "This is free software, and you are welcome to redistribute it\n"
	    << "under certain conditions. See the licence you obtained with\n"
	    << "the code\n"
	       "Git Checkout Hash " << GITHASH << "\n\n";
  ////////////////////////PROGRAM OPTIONS!!!!!!!!!!!!!!!!!!!!!!!
  try 
    {
      po::options_description allopts("General Options"), loadopts("Load Config File Options");
      
      allopts.add_options()
	("help", "Produces this message.")
	("out-config-file,o", 
	 po::value<string>()->default_value("config.out.xml.bz2"), 
	 "Configuration output file.")
	("random-seed,s", po::value<unsigned int>(),
	 "Random seed for generator.")
	("rescale-T,r", po::value<double>(), 
	 "Rescale kinetic temperature to this value.")
	("zero-momentum,Z", "Zero the momentum.")
	("zero-com", "Zero the centre of mass.")
	("zero-vel", po::value<size_t>(), "Set all particles velocity component (arg) to zero.")
	("set-com-vel", po::value<std::string>(), 
	 "Sets the velocity of the COM of the system (format x,y,z no spaces).")
	("mirror-system,M",po::value<unsigned int>(), 
	 "Mirror the particle co-ordinates and velocities. Argument is "
	 "dimension to reverse/mirror.")
	("round", "Output the XML config file with one less digit to aid in removing rounding errors for"
	 " test systems.")
	("uncompressed", "Output the XML config file without bzip compression.")
	;

      loadopts.add_options()
	("config-file", po::value<string>(), 
	 "Config file to initialise from (Non packer mode).")
	;
       
      allopts.add(loadopts);
      allopts.add(CIPPacker::getOptions());

      po::positional_options_description p;
      p.add("config-file", 1);
      
      po::variables_map vm;
      po::store(po::command_line_parser(argc, argv).
		options(allopts).positional(p).run(), vm);
      po::notify(vm);
      
      if ((vm.count("help") 
	   || (!vm.count("packer-mode") 
	       && !vm.count("config-file")))
	  && !vm.count("packer-mode-help"))
	{
	  cout << "Usage : dynamod <OPTION>...[CONFIG FILE]\n"
	       << " Modifies a config file (if provided) or generates a new"
	    " config file\n" << allopts << "\n";

	  return 1;
	}

      if (vm.count("random-seed"))
	sim.setRandSeed(vm["random-seed"].as<unsigned int>());
      
      ////////////////////////Simulation Initialisation!!!!!!!!!!!!!
      //Now load the config

      {
	std::string fileName(vm["out-config-file"].as<std::string>());
	if (vm.count("uncompressed") 
	    && (std::string(fileName.end()-4, fileName.end()) == ".bz2"))
	  M_throw() << "You should not use a .bz2 extension for uncompressed"
	    " comfig files";
      }

      if (!vm.count("config-file"))
	{
	  CIPPacker plug(vm, &sim);
	  plug.initialise();

	  std::cout << "\nMain: Finialising the packing routines";

	  //We don't zero momentum and rescale for certain packer modes
	  if ((vm["packer-mode"].as<size_t>() != 23)
	      && (vm["packer-mode"].as<size_t>() != 25))
	    {
	      CInputPlugin(&sim, "Rescaler").zeroMomentum();
	      CInputPlugin(&sim, "Rescaler").rescaleVels(1.0);
	    }
	  sim.configLoaded();
	}
      else
	sim.loadXMLfile(vm["config-file"].as<string>());
  
      sim.setTrajectoryLength(0);

      CIPPacker(vm, &sim).processOptions();

      sim.initialise();      
      
      //Here we modify the sim accordingly      

      if (vm.count("zero-momentum"))
	CInputPlugin(&sim, "MomentumZeroer")
	  .zeroMomentum();	

      if (vm.count("zero-com"))
	CInputPlugin(&sim, "CentreOfMassZeroer")
	  .zeroCentreOfMass();	

      if (vm.count("rescale-T"))
	CInputPlugin(&sim, "Rescaler")
	  .rescaleVels(vm["rescale-T"].as<double>());

      if (vm.count("mirror-system"))
	CInputPlugin(&sim, "Mirrorer").
	  mirrorDirection(vm["mirror-system"].as<unsigned int>());

      if (vm.count("set-com-vel"))
	{
	  boost::tokenizer<boost::char_separator<char> > 
	    tokens(vm["set-com-vel"].as<std::string>(), boost::char_separator<char>(","));
	  
	  boost::tokenizer<boost::char_separator<char> >::iterator details_iter = tokens.begin();

	  Vector vel(0,0,0);

	  if (details_iter == tokens.end()) M_throw() << "set-com-vel requires 3 components";
	  vel[0] = boost::lexical_cast<double>(*(details_iter++));
	  if (details_iter == tokens.end()) M_throw() << "set-com-vel requires 3 components";	  
	  vel[1] = boost::lexical_cast<double>(*(details_iter++));
	  if (details_iter == tokens.end()) M_throw() << "set-com-vel requires 3 components";
	  vel[2] = boost::lexical_cast<double>(*(details_iter));
	  
	  CInputPlugin(&sim, "velSetter")
	    .setCOMVelocity(vel);
	}

      if (vm.count("zero-vel"))
	CInputPlugin(&sim, "Vel-Component-Zeroer").
	  zeroVelComp(vm["zero-vel"].as<size_t>());


      //Write out now we've changed the system
      sim.getHistory() << "\nconfigmod run as so\n";
      for (int i = 0; i< argc; i++)
	sim.getHistory() << argv[i] << " ";
      sim.getHistory() << "\nGIT hash " << GITHASH;
      cout << "\nWriting out configuration";
      sim.writeXMLfile(vm["out-config-file"].as<string>(), vm.count("round"));
      cout << "\n";
    }
  catch (std::exception &cep)
    {
      std::cout << "\nReached Main Error Loop"
		<< "\nOutputting results so far and shutting down"
		<< "\nBad configuration written to config.error.xml"
		<< cep.what();
      
      try
	{ sim.writeXMLfile("config.error.xml.bz2"); }
      catch (std::exception &cep)
	{
	  std::cout << "\nFailed to output error config"
		    << cep.what() << "\n";
	}

      std::cout << "\n";

      return 1;
    }

  return 0;
}
