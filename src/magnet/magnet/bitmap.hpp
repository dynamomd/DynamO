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

#pragma once

#include <magnet/PNG.hpp>

namespace magnet {
  namespace image {
    namespace detail {
      template<class T>     
      inline void w(std::ostream& os, const T& obj)
      { os.write((const char*)&obj, sizeof(T)); }

      class bitmap_information_header
      {
      public:
	inline bitmap_information_header(size_t width, size_t height): 
	  _type(19778),
	  _fileSize(54 + ((width * 3 + 3) & 0xFFFC) * height),
	  _reserved1(0), _reserved2(0), _offset_bits(54),
	  _headerSize(40), _width(width), _height(height), _planes(1), _bitdepth(24),
	  _compression(0), _compressed_file_size(((width * 3 + 3) & 0xFFFC) * height),
	  _xres(2835), _yres(2835), _palleteSize(0), _importantColors(0) {}

	inline friend std::ostream& operator<<(std::ostream& os, 
					       const bitmap_information_header& b)
	{
	  w(os, b._type);
	  w(os, b._fileSize);
	  w(os, b._reserved1);
	  w(os, b._reserved2);
	  w(os, b._offset_bits);

	  w(os, b._headerSize);
	  w(os, b._width);
	  w(os, b._height);
	  w(os, b._planes); 
	  w(os, b._bitdepth);
	  w(os, b._compression); 
	  w(os, b._compressed_file_size);
	  w(os, b._xres);
	  w(os, b._yres);
	  w(os, b._palleteSize);
	  w(os, b._importantColors);

	  return os;
	}

      private:
	//File header - 14 bytes
	const uint16_t _type;
	const uint32_t _fileSize;
	const uint16_t _reserved1;
	const uint16_t _reserved2;
	const uint32_t _offset_bits;

	//The BMP header - 40 bytes
	const uint32_t _headerSize;
	const int32_t _width;
	const int32_t _height;
	const uint16_t _planes;
	const uint16_t _bitdepth;
	const uint32_t _compression;
	const uint32_t _compressed_file_size;
	const int32_t _xres;
	const int32_t _yres;
	const uint32_t _palleteSize;
	const uint32_t _importantColors;
      };
    }
    
    template<png::ColorType color>
    inline void writeBMPFile(const std::string& filename,
			     const std::vector<png::Pixel<color> >& image, 
			     size_t width, size_t height)
    {
      std::ofstream bmpFile(filename.c_str(), std::fstream::binary);
      bmpFile << detail::bitmap_information_header(width, height);
      
      size_t rowpadding = ((width * 3 + 3) & 0xFFFC) - width * 3;
      for (size_t y(0); y < height; ++y)
	{
	  for (size_t x(0); x < width; ++x)
	    {
	      detail::w(bmpFile, image[y * width + x].blue());
	      detail::w(bmpFile, image[y * width + x].green());
	      detail::w(bmpFile, image[y * width + x].red());
	    }
	  
	  for (size_t xpad(0); xpad < rowpadding; ++xpad)
	    detail::w(bmpFile,int8_t(0));
	}
    }
  }
}
