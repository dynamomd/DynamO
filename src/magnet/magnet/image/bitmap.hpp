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
#include <magnet/image/pixel.hpp>

namespace magnet {
  namespace image {
    namespace detail {
      /*! \brief Function to write the passed variable in binary to
       * the passed output stream.
       */
      template<class T>     
      inline void write(std::ostream& os, const T& obj)
      { os.write((const char*)&obj, sizeof(T)); }

      /*! \brief This class outputs the binary magic header for bitmap
       * files.
       */
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
	  write(os, b._type);
	  write(os, b._fileSize);
	  write(os, b._reserved1);
	  write(os, b._reserved2);
	  write(os, b._offset_bits);
	  write(os, b._headerSize);
	  write(os, b._width);
	  write(os, b._height);
	  write(os, b._planes); 
	  write(os, b._bitdepth);
	  write(os, b._compression); 
	  write(os, b._compressed_file_size);
	  write(os, b._xres);
	  write(os, b._yres);
	  write(os, b._palleteSize);
	  write(os, b._importantColors);

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
    
    /*! \brief Helper function to write out a collection of pixels as a bmp file.
     */
    template<ColorType color>
    inline void writeBMPFile(const std::string& filename,
			     const std::vector<Pixel<color> >& image, 
			     size_t width, size_t height)
    {
      std::ofstream bmpFile(filename.c_str(), std::fstream::binary);
      bmpFile << detail::bitmap_information_header(width, height);
      
      size_t rowpadding = ((width * 3 + 3) & 0xFFFC) - width * 3;
      for (size_t y(0); y < height; ++y)
	{
	  for (size_t x(0); x < width; ++x)
	    {
	      detail::write(bmpFile, image[y * width + x].blue());
	      detail::write(bmpFile, image[y * width + x].green());
	      detail::write(bmpFile, image[y * width + x].red());
	    }
	  
	  for (size_t xpad(0); xpad < rowpadding; ++xpad)
	    detail::write(bmpFile,int8_t(0));
	}
    }
  }
}
