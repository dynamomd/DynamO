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

#include <magnet/image/PNG.hpp>

namespace po = boost::program_options;

std::vector<magnet::image::Pixel<magnet::image::RGB> >  pixels;


int main(int argc, char *argv[])
{
  try {
    po::options_description allopts;
    allopts.add_options()
      ("help,h", "Produces this message.")
      ("font-file", po::value<std::string>(), "The font file to generate the texture from.")
      ("output-texture", po::value<std::string>()->default_value("out.png"), "The generated font texture.")
      ("size", po::value<size_t>()->default_value(48), "The font size in pixels.")
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
      throw std::runtime_error("Failed to initialise the freetype library.");

    FT_Face face;
    switch (FT_New_Face(ftlibrary, vm["font-file"].as<std::string>().c_str(), 0, &face))
      {
      case FT_Err_Unknown_File_Format:
	throw std::runtime_error("Failed to load the font file, unsupported file format.");
      default:
	throw std::runtime_error("Unknown Error: Failed to load the font file.");
      case 0: break;
      }

    if (!face)
      throw std::runtime_error("Face pointer is invalid.");

    if (!(face->face_flags & FT_FACE_FLAG_SCALABLE))
      throw std::runtime_error("The font is not a scalable font.");      
    
    if (FT_Set_Char_Size(face, 0, vm["size"].as<size_t>() * 64, 72, 72))
      throw std::runtime_error("Could not set font size");

    std::string input("String to render!");

    size_t lineheight = face->height;
    size_t max_charwidth = face->max_advance_width;
    size_t linewidth = max_charwidth * input.size();

    std::vector<magnet::image::Pixel<magnet::image::RGB> >  pixels;
    pixels.resize(linewidth * lineheight);

    size_t pen_x(0);
    size_t pen_y(0);

    for (std::string::const_iterator cPtr = input.begin(); 
	 cPtr != input.end(); ++cPtr)
      {
	if (FT_Load_Char(face, *cPtr, FT_LOAD_RENDER))
	  continue;

	size_t rows = face->glyph->bitmap.rows;
	size_t width = face->glyph->bitmap.width;

	for (size_t y(0); y < rows; ++y)
	  for (size_t x(0); x < width; ++x)
	    {
	      size_t xpos = pen_x + x + face->glyph->bitmap_left;
	      size_t ypos = pen_y + y + lineheight - face->glyph->bitmap_top;
	      if ((xpos < linewidth) && (ypos < lineheight))
		{
		  size_t value = face->glyph->bitmap.buffer[y * width + x];
		  pixels[ypos * linewidth + xpos]
		    = magnet::image::Pixel<magnet::image::RGB>(value,value,value);
		}
	    }
	pen_x += face->glyph->advance.x / 64;
      }

    magnet::image::writePNGFile(vm["output-texture"].as<std::string>(), pixels, 
				linewidth, lineheight);

  } catch (std::exception& err)
    {
      std::cout << "Exception caught : " 
		<< err.what()
		<< std::endl; 
      return 1;
    }
  
  return 0;
}
