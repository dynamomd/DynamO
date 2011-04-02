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

	std::vector<unsigned char> buffer(_width * _height * _depth);
	file.read(reinterpret_cast<char*>(&buffer[0]), buffer.size());
	
	//Now we must dilate the array, placing the value in the alpha
	//channel and calculate gradients
	std::vector<unsigned char> buffer2(4 * _width * _height * _depth);
	for (size_t i(0); i < buffer.size(); ++i)
	  {
	    buffer2[4 * i + 0] = buffer[i];
	    buffer2[4 * i + 1] = buffer[i];
	    buffer2[4 * i + 2] = buffer[i];
	    buffer2[4 * i + 3] = (buffer[i] > 16) ? 16 : 0;
	  }

	glTexSubImage3D(GL_TEXTURE_3D, 0, //level of detail
			0, 0, 0,//offset
			_width, _height, _depth, //size of the texture
			_format, _type, &buffer2[0]
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
      GLuint _handle;
      bool _valid;
      size_t _width, _height, _depth;
      GLenum _format;
      GLenum _type;
    };
  }
}
