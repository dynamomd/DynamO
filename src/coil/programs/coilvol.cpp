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
	("data-file", po::value<std::vector<string> >(), "Data file to display.")
	("x-elements,x", po::value<size_t>(), "Number of x volume elements")
	("y-elements,y", po::value<size_t>(), "Number of y volume elements")
	("z-elements,z", po::value<size_t>(), "Number of z volume elements")
	("data-size", po::value<size_t>()->default_value(1), 
	 "Size of each volume element (in bytes).")
	;

      po::positional_options_description p;
      p.add("data-file", -1);
      
      po::variables_map vm;
      po::store(po::command_line_parser(argc, argv).
		options(opts).positional(p).run(), vm);
      po::notify(vm);
      
      if (vm.count("help"))
	{
	  cout << "Usage : coilvol <OPTIONS>...[DATA FILE]\n"
	       << "Draws a raw volume data file using the coil library, "
	       << "you must set the data file, x, y, and z dimensions "
	       << "of the data set. Alternatively, just specify a list of TIFF files to be stacked." 
	       << opts;
	  return 1;
	}
      
      magnet::ArgShare::getInstance().setArgs(argc, argv);

      coil::CoilMaster::parallel = false;
      coil::CoilRegister coil;
      std::shared_ptr<coil::CLGLWindow> window(new coil::CLGLWindow("Coil Volume Renderer : ", 1.0));
      coil.getInstance().addWindow(window);
      
      if (vm.count("data-file")) {
	std::vector<std::string> files = vm["data-file"].as<std::vector<std::string> >();
	std::shared_ptr<coil::RVolume> voldata(new coil::RVolume("Volume data"));
	window->addRenderObj(voldata);
	
	if (files.size() == 1)
	  {
	    size_t datasize[3] = {vm["x-elements"].as<size_t>(), vm["y-elements"].as<size_t>(), vm["z-elements"].as<size_t>()};
	    window->getGLContext()->queueTask(std::bind(&coil::RVolume::loadRawFile, voldata.get(), files[0], std::array<size_t, 3>{{datasize[0], datasize[1], datasize[2]}}, vm["data-size"].as<size_t>()));
	  }
	else
	  {
#ifdef COIL_TIFFSUPPORT
	    std::cout << "Loading " << vm.count("data-file") << " datafiles" << std::endl;
	    window->getGLContext()->queueTask(std::bind(&coil::RVolume::loadTIFFFiles, voldata.get(), files));
#else
	    M_throw() << "Loading multiple images is only supported if TIFF support is built in";
#endif
	  }
      }

      auto& master = coil.getInstance();
      while (master.main_loop_iter()) {};
    }
  catch (std::exception &cep)
    {
      std::cerr << "\nException caught in main()\n"
		<< cep.what();
      return 1;
    }

  return 0;
}
