#ifndef PNG_HPP
#define PNG_HPP

// C++
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>

// C
#include <inttypes.h>

// PNG
#include <png.h>

class PNGStreamUtils {
public:
  static void read(png_structp png, png_bytep data, png_size_t length) {
    static_cast<std::istream*>(png_get_io_ptr(png))->read(
							  reinterpret_cast<char*>(data), length);
  }

  static void write(png_structp png, png_bytep data, png_size_t length) {
    static_cast<std::ostream*>(png_get_io_ptr(png))->write(
							   reinterpret_cast<char*>(data), length);
  }

  static void flush(png_structp png) {
    static_cast<std::ostream*>(png_get_io_ptr(png))->flush();
  }
};

class PNGPixel {
public:
  PNGPixel() : r_(0), g_(0), b_(0), a_(0) {
  }

  PNGPixel(uint32_t value) : r_(0), g_(0), b_(0), a_(0) {
    convert(value);
  }

  uint8_t red() const {
    return r_;
  }

  uint8_t& red() {
    return r_;
  }

  uint8_t green() const {
    return g_;
  }

  uint8_t& green() {
    return g_;
  }

  uint8_t blue() const {
    return b_;
  }

  uint8_t& blue() {
    return b_;
  }

  uint8_t alpha() const {
    return a_;
  }

  uint8_t& alpha() {
    return a_;
  }

  uint32_t convert() const {
    return (r_ << 24) + (g_ << 16) + (b_ << 8) + a_;
  }

  void convert(uint32_t value) {
    r_ = uint8_t((value >> 24) & 0xff);
    g_ = uint8_t((value >> 16) & 0xff);
    b_ = uint8_t((value >>  8) & 0xff);
    a_ = uint8_t(value & 0xff);
  }

private:
  uint8_t r_, g_, b_, a_;
};

inline std::ostream& operator<<(std::ostream& os, const PNGPixel& pixel) {
  os << "[" << size_t(pixel.red()) << ", " << size_t(pixel.green()) <<
    ", " << size_t(pixel.blue()) << ", " << size_t(pixel.alpha()) <<
    "] (" << pixel.convert() << ")";

  return os;
}

class PNGImage {
public:
  static void readFile(const std::string& filename,
		       std::vector<uint32_t>& image, size_t& width, size_t& height) {

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
		    PNGStreamUtils::read);

    png_set_sig_bytes(png, pngHeaderSize);

    png_read_info(png, pngInfo);

    if(png_get_color_type(png, pngInfo) != PNG_COLOR_TYPE_RGBA &&
       png_get_color_type(png, pngInfo) != PNG_COLOR_TYPE_RGB) {

      std::stringstream strm;
      strm << "unsupported color type in '" << filename << "'";
      throw std::runtime_error(strm.str().c_str());
    }

    // strip unneeded alpha channel
    if(png_get_color_type(png, pngInfo) & PNG_COLOR_MASK_ALPHA)	
      png_set_strip_alpha(png);

    // update png information before reading
    png_read_update_info(png, pngInfo);

    width = png_get_image_width(png, pngInfo);
    height = png_get_image_height(png, pngInfo);
    size_t channels = png_get_channels(png, pngInfo);
    //			size_t bitDepth = png_get_bit_depth(png, pngInfo);

    /*
      std::cerr << "reading png image of size " << width << "x" <<
      height << " pixels (" << bitDepth << " bit/pixel, " <<
      channels << " channels)..." << std::endl;
    */

    size_t bytesPerRow = png_get_rowbytes(png, pngInfo);
    png_bytep* pngRows = new png_bytep[height];
    png_bytep pngData = new png_byte[height * bytesPerRow];

    size_t offset = 0;
    for(size_t row = 0; row < height; ++row) {
      pngRows[row] = pngData + offset;
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

    // convert image to vector of uints
    image.resize(height * width);

    for(size_t y = 0; y < height; ++y)
      for(size_t x = 0; x < width; ++x) {
	PNGPixel value;
	value.red() = pngData[bytesPerRow * y + x * channels];

	value.green() = pngData[bytesPerRow * y + x * channels + 1];

	value.blue() = pngData[bytesPerRow * y + x * channels + 2];

	image[y * width + x] = value.convert();
      }

    // all went well, we can clean up now
    png_destroy_read_struct(&png, &pngInfo, &pngEndInfo);
    delete[] pngData;
    delete[] pngRows;
  }

  static void writeFile(const std::string& filename,
			const std::vector<uint32_t>& image, size_t& width,
			size_t& height, int compression_level = 6, bool gldata = false) {

    if(image.size() != width * height)
      throw std::runtime_error("invalid input vector in PNGImage::writeFile");

    std::ofstream pngFile(filename.c_str(), std::fstream::binary);

    if(!pngFile.is_open()) {
      std::stringstream strm;
      strm << "failed to open file '" << filename << "'";
      throw std::runtime_error(strm.str().c_str());
    }

    pngFile.exceptions(std::fstream::badbit | std::fstream::failbit);

    /*
      const size_t pngHeaderSize = 8;
      png_byte pngHeader[pngHeaderSize];

      pngFile.read(reinterpret_cast<char*>(pngHeader), pngHeaderSize);

      if(png_sig_cmp(pngHeader, 0, pngHeaderSize)) {
      std::stringstream strm;
      strm << "failed to read '" << filename << "': not a png file";
      throw std::runtime_error(strm.str().c_str());
      }
    */
    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING,
					      NULL, NULL, NULL);
			
    png_set_compression_level(png, compression_level);

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
		     PNGStreamUtils::write, PNGStreamUtils::flush);

    //			png_set_sig_bytes(png, pngHeaderSize);

    png_set_IHDR(png, pngInfo, width, height, 8, PNG_COLOR_TYPE_RGB,
		 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE,
		 PNG_FILTER_TYPE_BASE);

    png_write_info(png, pngInfo);

    size_t channels = 3;		// RGB only
    size_t rowSize = channels * width;
    png_bytep* pngRows = new png_bytep[height];
    png_bytep pngData = new png_byte[height * rowSize];

    size_t offset = 0;
    for(size_t row = 0; row < height; ++row) {
      pngRows[row] = pngData + offset;
      offset += rowSize;
    }

    // convert vector of uints to image
    if (!gldata)
      for(size_t y = 0; y < height; ++y)
	for(size_t x = 0; x < width; ++x) 
	  {
	    PNGPixel value;
	    value.convert(image[y * width + x]);
	    pngData[rowSize * y + x * channels] = value.red();
	    pngData[rowSize * y + x * channels + 1] = value.green();
	    pngData[rowSize * y + x * channels + 2] = value.blue();
	  }
    else
      for(size_t y = 0; y < height; ++y)
	for(size_t x = 0; x < width; ++x) 
	  {
	    uint32_t pixel = image[(height - 1 - y) * width + x];
	    pngData[rowSize * y + x * channels] = pixel & 0xff;
	    pngData[rowSize * y + x * channels + 1] = (pixel >> 8) & 0xff;
	    pngData[rowSize * y + x * channels + 2] = (pixel >> 16) & 0xff;
	  }


    // set up error handling the hard way
    if(setjmp(png_jmpbuf(png))) {
      png_destroy_write_struct(&png, &pngInfo);
      delete[] pngData;
      delete[] pngRows;
      throw std::runtime_error("libpng: failed to write image");
    }

    png_write_image(png, pngRows);
    png_write_end(png, NULL);

    // all went well, we can clean up now
    png_destroy_write_struct(&png, &pngInfo);
    delete[] pngData;
    delete[] pngRows;
  }
};

#endif // PNG_HPP
