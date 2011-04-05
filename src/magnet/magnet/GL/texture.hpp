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
      GLint _internalFormat;
      const GLenum _texType;
    };

    class Texture1D: public TextureBasic
    {
    public:
      Texture1D():TextureBasic(GL_TEXTURE_1D) {}
      
      inline void init(size_t width, GLint internalformat = GL_RGBA8)
      {
	_width = width; _internalFormat = internalformat;
	TextureBasic::init();
	bind(0);
	parameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	parameter(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);

	glTexImage1D(_texType, 0, _internalFormat, 
		     _width,
		     //Border is off
		     0,
		     //The following values are not used as no data is
		     //passed here, use subImage for that.
		     GL_RGB, GL_UNSIGNED_BYTE, NULL); 
      }

      inline void subImage(const std::vector<GLubyte>& data, GLenum pixelformat,
			   GLint xoffset = 0, GLint width = -1, GLint level = 0)
      { 
	//If there are negative values, assume they mean to use the
	//full space
	if (width  < 0) width  = _width;
	if (xoffset < 0) M_throw() << "x offset is negative";
	if (xoffset + width > _width) M_throw() << "Texture write x overrun";

	bind(0);
	glTexSubImage1D(_texType, level, xoffset, width,
			pixelformat, GL_UNSIGNED_BYTE, &data[0]);
      }

      inline const GLint& getWidth() const { return _width; }
    private:
      
      GLint _width;
    };


    class Texture3D: public TextureBasic
    {
    public:
      Texture3D():TextureBasic(GL_TEXTURE_3D) {}
      
      inline void init(size_t width, size_t height, size_t depth, 
		       GLint internalformat = GL_RGBA8)
      {
	_width = width; _height = height; _depth = depth;
	_internalFormat = internalformat;
	TextureBasic::init();
	bind(0);

	parameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	parameter(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	parameter(GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

	glTexImage3D(_texType, 0, _internalFormat, 
		     _width, _height, _depth, 
		     //Border is off
		     0,
		     //The following values are not used as no data is
		     //passed here, use subImage for that.
		     GL_RGB, GL_UNSIGNED_BYTE, NULL); 
      }

      inline void subImage(const std::vector<GLubyte>& data, GLenum pixelformat,
			   GLint xoffset = 0, GLint yoffset = 0, GLint zoffset = 0,
			   GLint width = -1, GLint height = -1, GLint depth = -1,
			   GLint level = 0)
      { 
	//If there are negative values, assume they mean to use the
	//full space
	if (width  < 0) width  = _width;
	if (height < 0) height = _height;
	if (depth  < 0) depth  = _depth;

	if (xoffset < 0) M_throw() << "x offset is negative";
	if (yoffset < 0) M_throw() << "y offset is negative";
	if (zoffset < 0) M_throw() << "z offset is negative";
	
	if (xoffset + width > _width) M_throw() << "Texture write x overrun";
	if (yoffset + height > _height) M_throw() << "Texture write y overrun";
	if (zoffset + depth > _depth) M_throw() << "Texture write z overrun";

	bind(0);
	glTexSubImage3D(_texType, level, xoffset, yoffset, zoffset, width, height, depth,
			pixelformat, GL_UNSIGNED_BYTE, &data[0]);
      }

      inline const GLint& getWidth() const { return _width; }
      inline const GLint& getHeight() const { return _height; }
      inline const GLint& getDepth() const { return _depth; }

    private:
      
      GLint _width, _height, _depth;
    };

  }
}
