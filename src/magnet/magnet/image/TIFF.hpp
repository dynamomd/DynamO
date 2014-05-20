/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
    Copyright (C) 2013  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#pragma once
#include <tiffio.h>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <algorithm>

namespace magnet {
  namespace image {
    /*Pixels have four 8 bit channels, but may also be viewed as a single
      32 bit int. Here we use the dreaded union and bit field to make it
      appear as either.*/
    struct Pixel {
      uint8 r;
      uint8 g;
      uint8 b;
      uint8 a;
    };

    struct Image {
      uint32 width;
      uint32 height;
      std::vector<Pixel> pixels;
    };

    /* Function to load a single TIFF file as an Image. */
    Image loadTIFF(const std::string filename)
    {
      TIFF* tif = TIFFOpen(filename.c_str(), "r");
      if (!tif)
	throw std::runtime_error("Failed to open image file");
      Image img;
      TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &(img.width));
      TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &(img.height));
      img.pixels.resize(img.width * img.height);
      if (!TIFFReadRGBAImage(tif, img.width, img.height, reinterpret_cast<uint32*>(&img.pixels[0]), 0))
	throw std::runtime_error("Failed to read image data");
      TIFFClose(tif);
      return img;
    }

    struct Volume {
      Volume():width(0), height(0), depth(0) {}
      size_t width;
      size_t height;
      size_t depth;
      std::vector<Pixel> pixels;
    };


    /* Function to load a list of TIFF files, and stack them into a Volume image. */
    Volume loadTIFFStack(std::vector<std::string> filenames)
    {
      //Sort the filenames, as this order determines the stack
      std::sort(filenames.begin(), filenames.end());
  
      Volume vol;
      for (std::vector<std::string>::iterator fptr = filenames.begin(); fptr != filenames.end(); ++fptr)
	{
	  Image img = loadTIFF(*fptr);
	  if (vol.width == 0) {
	    vol.width = img.width;
	    vol.height = img.height;
	  }
      
	  if ((img.width != vol.width) || (img.height != vol.height))
	    throw std::runtime_error("Images have varying dimensions");
	  const size_t last_size = vol.pixels.size();
	  vol.pixels.resize(vol.pixels.size() + vol.width * vol.height);
	  std::copy(img.pixels.begin(), img.pixels.end(), vol.pixels.begin()+last_size);
	  ++vol.depth;
	}
      return vol;
    }
  }
}
