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
#include <iostream>
#include <stdexcept>
#include <boost/program_options.hpp>

#include <ft2build.h>
#include FT_FREETYPE_H

namespace po = boost::program_options;

int main(int argc, char *argv[])
{
  po::options_description allopts;
  allopts.add_options()
    ("help,h", "Produces this message.")
    ("font-file", po::value<std::string>(), "The font file to generate the texture from.")
    ("output-texture", po::value<std::string>()->default_value("out.raw"), "The generated font texture.")
    ;

  po::positional_options_description p;
  p.add("font-file", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).
	    options(allopts).positional(p).run(), vm);
  po::notify(vm);

  if (!vm.count("font-file") || vm.count("help"))
    {
      std::cout << "Usage : " << argv[0] << " <OPTIONS> [Font-File] \n"
		<< allopts;
      return 1;
    }

  FT_Library ftlibrary;
  if (FT_Init_FreeType(&ftlibrary))
    throw std::runtime_error("Failed to initialise the freetype library");

  FT_Face face;
  switch (FT_New_Face(ftlibrary, vm["font-file"].as<std::string>().c_str(), 0, &face))
    {
    case FT_Err_Unknown_File_Format:
      throw std::runtime_error("Failed to load the font file, unsupported file format!");
    default:
      throw std::runtime_error("Unknown Error: Failed to load the font file");
    case 0: break;
    }
}
