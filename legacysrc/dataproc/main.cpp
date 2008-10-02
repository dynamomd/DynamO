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

#include <gtkmm/main.h>
#include <iostream>
#include "window.hpp"
#include "../base/is_exception.hpp"
#include <boost/program_options.hpp>
using namespace std;
using namespace boost;
namespace po = boost::program_options;

int main(int argc, char *argv[])
{
  po::options_description allopts("");
  
  allopts.add_options() 	
    ("data-file", po::value<string>(), "Data file to initialise from")
    ;
  
  po::positional_options_description p;
  p.add("data-file", 1);
  
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).
	    options(allopts).positional(p).run(), vm);
  po::notify(vm);

  try {
    Gtk::Main kit(argc, argv);
    

    if (vm.count("data-file"))   
      {
	MainWindow window(vm["data-file"].as<std::string>());
	//Shows the window and returns when it is closed.
	Gtk::Main::run(window);	
      }
    else
      {
	MainWindow window;
	//Shows the window and returns when it is closed.
	Gtk::Main::run(window);
      }
    
    return 0;
  } catch (DYNAMO::Exception& cep)
    {
      std::cerr << cep.what() << "\n";
    }
}
