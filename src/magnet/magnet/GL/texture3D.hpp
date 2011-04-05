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

    class TextureBasic
    {
    public:
      inline TextureBasic(GLenum texType):
	_valid(false),
	_texType(texType)
      {}

      inline ~TextureBasic() { deinit(); }

      inline void init()
      {
	if (_valid) M_throw() << "Already init()ed!";
	glGenTextures(1, &_handle);
	_valid = true;
      }

      inline void deinit()
      {
	if (_valid)
	  {
	    glDeleteTextures(1, &_handle);
	    _valid = false;
	  }
      }

      inline void bind(int unit)
      {
	glActiveTextureARB(GL_TEXTURE0 + unit);
	glBindTexture(_texType, _handle);
      }

      inline void parameter(GLenum paramname, GLint param)
      { bind(0); glTexParameteri(_texType, paramname, param); }

      inline void parameter(GLenum paramname, GLfloat param)
      { bind(0); glTexParameterf(_texType, paramname, param); }

    protected:
      
      GLuint _handle;
      bool _valid;
      GLenum _format;
      GLenum _type;
      GLint _internalFormat;
      const GLenum _texType;
    };

    class Texture3D: public TextureBasic
    {
    public:
      Texture3D():TextureBasic(GL_TEXTURE_3D) {}
      
      inline void init(size_t width, size_t height, size_t depth, 
		       GLint internalformat = GL_RGBA8, 
		       GLenum format = GL_RGBA, GLenum type = GL_UNSIGNED_BYTE)
      {
	_width = width; _height = height; _depth = depth;
	_format = format; _type = type; _internalFormat = internalformat;
	TextureBasic::init();
	bind(0);

	parameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	parameter(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	parameter(GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
	image();
      }

      inline void image(GLint level=0, const GLvoid* data = NULL, GLint border = 0)
      { 
	bind(0); 
	glTexImage3D(_texType, level, _internalFormat, 
		     _width, _height, _depth, 
		     border, _format, _type, data); 
      }

      inline void subImage(const GLvoid* data, 
			   GLint xoffset = 0, GLint yoffset = 0, GLint zoffset = 0,
			   GLint width = -1, GLint height = -1, GLint depth = -1,
			   GLint level = 0)
      { 
	if (width  < 0) width  = _width;
	if (height < 0) height = _height;
	if (depth  < 0) depth  = _depth;
	
	if (xoffset + width > _width) M_throw() << "Texture write x overrun";
	if (yoffset + height > _height) M_throw() << "Texture write y overrun";
	if (zoffset + depth > _depth) M_throw() << "Texture write z overrun";

	bind(0);
	glTexSubImage3D(_texType, level, xoffset, yoffset, zoffset, width, height, depth,
			_format, _type, data);
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

	subImage(&buffer[0]);
      }

    private:
      //This stuff is for volume rendering
      size_t coordCalc(int x, int y, int z) const
      { 
	x = std::min(std::max(x, 0), _width  - 1); 
	y = std::min(std::max(y, 0), _height - 1); 
	z = std::min(std::max(z, 0), _depth  - 1); 
	return x + _width * (y + _height * z);
      }
      std::vector<unsigned char> calcVolData(const std::vector<unsigned char>& buff)
      {
	std::vector<unsigned char> unsmoothed(4 * _width * _height * _depth);

	for (GLint z(0); z < _depth; ++z)
	  for (GLint y(0); y < _height; ++y)
	    for (GLint x(0); x < _width; ++x)
	      {
		Vector sample1(buff[coordCalc(x - 2, y, z)],
			       buff[coordCalc(x, y - 2, z)],
			       buff[coordCalc(x, y, z - 2)]);

		Vector sample2(buff[coordCalc(x + 2, y, z)],
			       buff[coordCalc(x, y + 2, z)],
			       buff[coordCalc(x, y, z + 2)]);
		
		//Note, we store the negative gradient (we point down
		//the slope)
		Vector grad = sample1 - sample2;

		float nrm = grad.nrm();
		if (nrm > 0) grad /= nrm;
		
		size_t coord = x + _width * (y + _height * z);
		unsmoothed[4 * coord + 0] = uint8_t((grad[0] * 0.5 + 0.5) * 255);
		unsmoothed[4 * coord + 1] = uint8_t((grad[1] * 0.5 + 0.5) * 255);
		unsmoothed[4 * coord + 2] = uint8_t((grad[2] * 0.5 + 0.5) * 255);
		unsmoothed[4 * coord + 3] = buff[coordCalc(x, y, z)];
	      }

	return unsmoothed;
      }
      
      GLint _width, _height, _depth;
    };
  }
}
