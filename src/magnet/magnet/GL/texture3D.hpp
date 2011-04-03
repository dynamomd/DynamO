/*    DYNAMO:- Event driven molecular dynamics simulator 
 *    http://www.marcusbannerman.co.uk/dynamo
 *    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>
 *
 *    This program is free software: you can redistribute it and/or
 *    modify it under the terms of the GNU General Public License
 *    version 3 as published by the Free Software Foundation.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include <magnet/GL/detail/filter.hpp>
#include <magnet/exception.hpp>
#include <fstream>
#include <algorithm>

namespace magnet {
  namespace GL {

    class Texture3D
    {
    public:
      Texture3D():_valid(false) {}
      ~Texture3D() { deinit(); }

      inline void init(size_t width, size_t height, size_t depth, 
		       GLint internalformat = GL_RGBA8, 
		       GLenum format = GL_RGBA, GLenum type = GL_UNSIGNED_BYTE)
      {
	_width = width; _height = height; _depth = depth;
	_format = format; _type = type;
	if (_valid) M_throw() << "Already init()ed!";

	glGenTextures(1, &_handle);
	glBindTexture(GL_TEXTURE_3D, _handle);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
	glTexImage3D(GL_TEXTURE_3D, 0, internalformat, _width, _height, _depth, 0, 
		     format, type, NULL);

	_valid = false;
      }

      inline void bind(int unit)
      {
	glActiveTextureARB(GL_TEXTURE0 + unit);
	glBindTexture(GL_TEXTURE_3D, _handle);
      }

      inline void readFromRawFile(std::string filename)
      {
	std::ifstream file(filename.c_str(), std::ifstream::binary);
	
	std::vector<unsigned char> inbuffer(_width * _height * _depth);
	file.read(reinterpret_cast<char*>(&inbuffer[0]), inbuffer.size());
	
	//Sphere test pattern
	//for (size_t z(0); z < _depth; ++z)
	//  for (size_t y(0); y < _height; ++y)
	//    for (size_t x(0); x < _width; ++x)
	//      inbuffer[x + _width * (y + _height * z)] 
	//	= (std::sqrt(  std::pow(x - _width / 2.0, 2) 
	//		     + std::pow(y - _height / 2.0, 2) 
	//		     + std::pow(z - _depth / 2.0, 2))
	//		     < 122.0) ? 255.0 : 0;

	std::vector<unsigned char> buffer = calcVolData(inbuffer);

	glTexSubImage3D(GL_TEXTURE_3D, 0, //level of detail
			0, 0, 0,//offset
			_width, _height, _depth, //size of the texture
			_format, _type, &buffer[0]
			);
      }
      
      inline void deinit()
      {
	if (_valid)
	  {
	    glDeleteTextures(1, &_handle);
	    _valid = false;
	  }
      }

    private:
      //This stuff is for volume rendering
      class ClampedVector
      {
      public:
	ClampedVector(int w, int h, int d):
	  _width(w), _height(h), _depth(d)
	{}
	
	size_t operator()(int x, int y, int z) const
	{
	  size_t coord = std::min(_width - 1, std::max(x, 0))
	    + _width * (std::min(_height - 1, std::max(y, 0))
			+ _height * std::min(_depth - 1, std::max(z, 0)));
	  return coord;
	}

      private:
	int _width, _height, _depth;
      };

      std::vector<unsigned char> calcVolData(const std::vector<unsigned char>& buff)
      {
	std::vector<unsigned char> unsmoothed(4 * _width * _height * _depth);
	ClampedVector coordc(_width, _height, _depth);

	for (size_t z(0); z < _depth; ++z)
	  for (size_t y(0); y < _height; ++y)
	    for (size_t x(0); x < _width; ++x)
	      {
		Vector sample1(buff[coordc(x - 2, y, z)],
			       buff[coordc(x, y - 2, z)],
			       buff[coordc(x, y, z - 2)]);

		Vector sample2(buff[coordc(x + 2, y, z)],
			       buff[coordc(x, y + 2, z)],
			       buff[coordc(x, y, z + 2)]);
		
		//Note, we store the negative gradient (we point down
		//the slope)
		Vector grad = sample1 - sample2;

		float nrm = grad.nrm();
		if (nrm > 0) grad /= nrm;
		
		size_t coord = x + _width * (y + _height * z);
		unsmoothed[4 * coord + 0] = uint8_t((grad[0] * 0.5 + 0.5) * 255);
		unsmoothed[4 * coord + 1] = uint8_t((grad[1] * 0.5 + 0.5) * 255);
		unsmoothed[4 * coord + 2] = uint8_t((grad[2] * 0.5 + 0.5) * 255);
		unsmoothed[4 * coord + 3] = buff[coordc(x, y, z)];
	      }

//	std::vector<unsigned char> retval(4 * _width * _height * _depth);
//
//	int n = 5;
//	float count = n*n*n;
//	n = (n-1)/2;
//	//Now smooth the gradients
//	for (size_t z(0); z < _depth; ++z)
//	  for (size_t y(0); y < _height; ++y)
//	    for (size_t x(0); x < _width; ++x)
//	      {
//		size_t coord = x + _width * (y + _height * z);
//		//Store the actual val
//		retval[4 * coord + 3] = unsmoothed[4 * coord + 3];
//
//		//Now smooth the gradients
//		for (size_t c(0); c < 3; ++c)
//		  {
//		    float sum = 0;
//		    for (int i(-n); i < n; ++i)
//		      for (int j(-n); j < n; ++j)
//			for (int k(-n); k < n; ++k)
//			  sum += unsmoothed[4 * coordc(x+i, y+j, z+k) + c];
//		    retval[4 * coord + c] = uint8_t(sum * 255.0 / count);
//		  }
//	      }
//
	return unsmoothed;
      }
      
      GLuint _handle;
      bool _valid;
      size_t _width, _height, _depth;
      GLenum _format;
      GLenum _type;
    };
  }
}
