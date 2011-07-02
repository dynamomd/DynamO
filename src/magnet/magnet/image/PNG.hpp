/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>
    Copyright (C) 2011  Severin Strobl <-->

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
// C++
#include <fstream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>

// PNG
#define PNG_SKIP_SETJMP_CHECK
#include <png.h>

#include <magnet/image/pixel.hpp>

namespace magnet {
  namespace image {
    namespace detail {
      //! \brief Stream read function to pass as parameter to the PNG library
      inline void read(png_structp png, png_bytep data, png_size_t length) 
      { static_cast<std::istream*>(png_get_io_ptr(png))->read(reinterpret_cast<char*>(data), length); }

      //! \brief Stream write function to pass as parameter to the PNG library
      inline void write(png_structp png, png_bytep data, png_size_t length) 
      { static_cast<std::ostream*>(png_get_io_ptr(png))->write(reinterpret_cast<char*>(data), length); }

      //! \brief Stream flush function to pass as parameter to the PNG library
      inline void flush(png_structp png)
      { static_cast<std::ostream*>(png_get_io_ptr(png))->flush(); }
    }

    /*! \brief Reads a PNG file and places the pixel data in the
     * passed array.
     *
     * \param filename The name of the file to read.
     * \param image The array of \ref Pixel data to read the image into.
     * \param width Variable used to return the width of the image.
     * \param height Variable used to return the height of the image.
     */
    template<ColorType CT>
    inline void readPNGFile(const std::string& filename,
			    std::vector<Pixel<CT> >& image, size_t& width,
			    size_t& height) {

      std::ifstream pngFile(filename.c_str(), std::fstream::binary);

      if(!pngFile.is_open()) {
	std::stringstream strm;
	strm << "failed to open file '" << filename << "'";
	throw std::runtime_error(strm.str().c_str());
      }

      pngFile.exceptions(std::fstream::badbit | std::fstream::failbit);

      const size_t pngHeaderSize = 8;
      png_byte pngHeader[pngHeaderSize];

      pngFile.read(reinterpret_cast<char*>(pngHeader), pngHeaderSize);

      if(png_sig_cmp(pngHeader, 0, pngHeaderSize)) {
	std::stringstream strm;
	strm << "failed to read '" << filename << "': not a png file";
	throw std::runtime_error(strm.str().c_str());
      }

      png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING,
					       NULL, NULL, NULL);

      if(!png)
	throw std::runtime_error("failed to allocate png_struct");

      png_infop pngInfo = png_create_info_struct(png);

      if(!pngInfo) {
	png_destroy_read_struct(&png, NULL, NULL);
	throw std::runtime_error("failed to allocate png_info_struct");
      }

      png_infop pngEndInfo = png_create_info_struct(png);

      if(!pngEndInfo) {
	png_destroy_read_struct(&png, &pngInfo, NULL);
	throw std::runtime_error("failed to allocate png_info_struct");
      }

      // set up error handling the hard way
      if(setjmp(png_jmpbuf(png))) {
	png_destroy_read_struct(&png, &pngInfo, &pngEndInfo);
	throw std::runtime_error("libpng: failed to set up io");
      }

      png_set_read_fn(png, static_cast<void*>(&pngFile),
		      detail::read);

      png_set_sig_bytes(png, pngHeaderSize);

      png_read_info(png, pngInfo);

      if(png_get_color_type(png, pngInfo) != PNG_COLOR_TYPE_RGBA &&
	 png_get_color_type(png, pngInfo) != PNG_COLOR_TYPE_RGB) {

	png_destroy_read_struct(&png, &pngInfo, &pngEndInfo);

	std::stringstream strm;
	strm << "unsupported color type in '" << filename << "'";
	throw std::runtime_error(strm.str().c_str());
      }

      bool needCopy = false;

      if(CT == RGBA && png_get_color_type(png, pngInfo) ==
	 PNG_COLOR_TYPE_RGB)
	needCopy = true;

      if(CT != RGBA && png_get_color_type(png, pngInfo) ==
	 PNG_COLOR_TYPE_RGBA) {

	std::stringstream strm;
	strm << "failed to read file '" << filename << "': " <<
	  "images uses alpha channel, parameter 'image' has to be " <<
	  " of type std::vector<magnet::image::Pixel<magnet::image::RGBA> >";
	throw std::runtime_error(strm.str().c_str());
      }

      width = png_get_image_width(png, pngInfo);
      height = png_get_image_height(png, pngInfo);
      size_t channels = png_get_channels(png, pngInfo);
      size_t bitDepth = png_get_bit_depth(png, pngInfo);

      if(channels < 3) {
	png_destroy_read_struct(&png, &pngInfo, &pngEndInfo);

	std::stringstream strm;
	strm << "failed to read '" << filename << "': unsupported " <<
	  "number of channels: " << channels;
	throw std::runtime_error(strm.str().c_str());
      }

      if(bitDepth != 8) {
	png_destroy_read_struct(&png, &pngInfo, &pngEndInfo);

	std::stringstream strm;
	strm << "failed to read '" << filename << "': invalid bit " <<
	  "depth: " << bitDepth;
	throw std::runtime_error(strm.str().c_str());
      }

#ifdef DEBUG
      std::cerr << "reading png image of size " << width << "x" <<
	height << " pixels (" << bitDepth << " bit/pixel, " <<
	channels << " channels)..." << std::endl;
#endif

      size_t bytesPerRow = png_get_rowbytes(png, pngInfo);
      png_bytep* pngRows = new png_bytep[height];
      png_bytep pngData = 0;

      if(needCopy)
	pngData = new png_byte[height * bytesPerRow];

      image.resize(height * width);

      size_t offset = 0;
      png_bytep basePointer = pngData;

      if(!needCopy)
	basePointer = reinterpret_cast<png_bytep>(&image[0]);

      for(size_t row = 0; row < height; ++row) {
	pngRows[row] = basePointer + offset;
	offset += bytesPerRow;
      }

      // set up error handling the hard way
      if(setjmp(png_jmpbuf(png))) {
	png_destroy_read_struct(&png, &pngInfo, &pngEndInfo);
	delete[] pngData;
	delete[] pngRows;
	throw std::runtime_error("libpng: failed to read image");
      }

      png_read_image(png, pngRows);
      png_read_end(png, pngEndInfo);

      if(needCopy) {
	for(size_t y = 0; y < height; ++y) {
	  for(size_t x = 0; x < width; ++x) {
	    Pixel<CT> value;

	    value.convert(pngData + bytesPerRow * y + x * channels);

	    image[y * width + x] = value;
	  }
	}
      }

      // all went well, we can clean up now
      png_destroy_read_struct(&png, &pngInfo, &pngEndInfo);
      delete[] pngData;
      delete[] pngRows;
    }

    /*! \brief Writes a PNG file using the pixel data in the passed
     * array.
     *
     * \param filename The name of the file to create/overwrite.
     * \param image The array of \ref Pixel data to write out.
     * \param width The width of the image.
     * \param height The height of the image.
     * \param compressionLevel The level of compression requested from the PNG library.
     * \param disableFiltering Prevent the PNG library from filtering the output.
     * \param flip Flips the vertical ordering of the image (to allow easy saving of OpenGL renderings).
     */
    template<ColorType CT>
    inline void writePNGFile(const std::string& filename,
			     std::vector<Pixel<CT> >& image, size_t width, size_t height,
			     int compressionLevel = PNG_COMPRESSION_TYPE_DEFAULT,
			     bool disableFiltering = false, bool flip = false) {

      if(image.size() != width * height) {
	std::stringstream strm;
	strm << "invalid input to writePNGFile(): " <<
	  "size mismatch of input vector (is " << image.size() <<
	  ", should be " << width << "x" << height << " = " <<
	  width * height;
	throw std::runtime_error(strm.str().c_str());
      }

      if(compressionLevel < 0 || compressionLevel > 9) {
	std::stringstream strm;
	strm << "invalid input to writePNGFile(): " <<
	  "valid compression levels range from 0 to 9 (default: " <<
	  PNG_COMPRESSION_TYPE_DEFAULT << ")";
	throw std::runtime_error(strm.str().c_str());
      }

      std::ofstream pngFile(filename.c_str(), std::fstream::binary);

      if(!pngFile.is_open()) {
	std::stringstream strm;
	strm << "failed to open file '" << filename << "'";
	throw std::runtime_error(strm.str().c_str());
      }

      pngFile.exceptions(std::fstream::badbit | std::fstream::failbit);

      png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING,
						NULL, NULL, NULL);

      if(!png)
	throw std::runtime_error("failed to allocate png_struct");

      png_infop pngInfo = png_create_info_struct(png);

      if(!pngInfo) {
	png_destroy_write_struct(&png, NULL);
	throw std::runtime_error("failed to allocate png_info_struct");
      }

      // set up error handling the hard way
      if(setjmp(png_jmpbuf(png))) {
	png_destroy_write_struct(&png, &pngInfo);
	throw std::runtime_error("libpng: failed to set up io");
      }

      png_set_write_fn(png, static_cast<void*>(&pngFile),
		       detail::write, detail::flush);

      png_set_IHDR(png, pngInfo, width, height, 8, CT,
		   PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
		   PNG_FILTER_TYPE_BASE);

      if(disableFiltering)
	png_set_filter(png, PNG_FILTER_TYPE_BASE, PNG_FILTER_NONE);

      if(compressionLevel != PNG_COMPRESSION_TYPE_DEFAULT)
	png_set_compression_level(png, compressionLevel);

      png_write_info(png, pngInfo);

      size_t bytesPerRow = png_get_rowbytes(png, pngInfo);
      png_bytep* pngRows = new png_bytep[height];

      if(image.size() * sizeof(Pixel<CT>) != (height * bytesPerRow)) 
	{
	  std::stringstream strm;
	  strm << "writePNGFile(): invalid size of input data";
	  throw std::runtime_error(strm.str().c_str());
	}

      size_t offset = 0;
      png_bytep basePointer = reinterpret_cast<png_bytep>(&image[0]);
			
      if (flip)
	for(int row = height - 1; row >= 0; --row) 
	  {
	    pngRows[row] = basePointer + offset;
	    offset += bytesPerRow;
	  }
      else
	for(size_t row = 0; row < height; ++row) 
	  {
	    pngRows[row] = basePointer + offset;
	    offset += bytesPerRow;
	  }

      // set up error handling the hard way
      if(setjmp(png_jmpbuf(png))) {
	png_destroy_write_struct(&png, &pngInfo);
	delete[] pngRows;
	throw std::runtime_error("libpng: failed to write image");
      }

      png_write_image(png, pngRows);
      png_write_end(png, NULL);

      // all went well, we can clean up now
      png_destroy_write_struct(&png, &pngInfo);
      delete[] pngRows;
    }
  }
}
