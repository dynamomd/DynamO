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

#include <coil/coilMaster.hpp>
#include <coil/clWindow.hpp>
#include <coil/RenderObj/Volume.hpp>
#include <magnet/arg_share.hpp>
#include <boost/program_options.hpp>
#include <iostream>

using namespace std;
using namespace boost;
namespace po = boost::program_options;

int main(int argc, char *argv[])
{
  try 
    {
      po::options_description opts("General Options");

      opts.add_options()
	("help,h", "Produces this message.")
	("data-file", po::value<string>(), "Data file to display.")
	("x-elements,x", po::value<size_t>(), "Number of x volume elements")
	("y-elements,y", po::value<size_t>(), "Number of y volume elements")
	("z-elements,z", po::value<size_t>(), "Number of z volume elements")
	("data-size", po::value<size_t>()->default_value(1), 
	 "Size of each volume element (in bytes).")
	;

      po::positional_options_description p;
      p.add("data-file", 1);
      
      po::variables_map vm;
      po::store(po::command_line_parser(argc, argv).
		options(opts).positional(p).run(), vm);
      po::notify(vm);
      
      if (vm.count("help") 
	  || !vm.count("data-file")
	  || !vm.count("x-elements")
	  || !vm.count("y-elements")
	  || !vm.count("z-elements"))
	{
	  cout << "Usage : coilvol <OPTIONS>...[DATA FILE]\n"
	       << "Draws a raw volume data file using the coil library, "
	       << "you must set the data file, x, y, and z dimensions "
	       << "of the data set." 
	       << opts;
	  return 1;
	}
      
      
      magnet::ArgShare::getInstance().setArgs(argc, argv);
      
      coil::CoilRegister coil;
            
      std::tr1::shared_ptr<coil::RVolume> 
	voldata(new coil::RVolume(vm["data-file"].as<std::string>()));

      std::tr1::shared_ptr<coil::CLGLWindow> 
	window(new coil::CLGLWindow("Coil Volume Renderer : ", 1.0));

      window->addRenderObj(voldata);
      
      coil.getInstance().addWindow(window);

      voldata->loadRawFile(vm["data-file"].as<std::string>(), 
			   vm["x-elements"].as<size_t>(),
			   vm["y-elements"].as<size_t>(),
			   vm["z-elements"].as<size_t>(),
			   vm["data-size"].as<size_t>());
      
      while (true) { window->simupdateTick(0); }
    }
  catch (std::exception &cep)
    {
      std::cerr << "\nException caught in main()\n"
		<< cep.what();
      return 1;
    }

  return 0;
}
