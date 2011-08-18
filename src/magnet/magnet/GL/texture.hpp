/*    dynamo:- Event driven molecular dynamics simulator 
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

#include <magnet/GL/context.hpp>
#include <magnet/exception.hpp>

namespace magnet {
  namespace GL {
    namespace detail {

      /*! \brief Generic interface for Texture objects. */
      class TextureBasic
      {
      protected:
	/*! \brief Constructor which requires the texture type. */
	inline TextureBasic(GLenum texType):
	  _valid(false),
	  _texType(texType)
	{}

	inline ~TextureBasic() { deinit(); }

	/*! \brief Allocates the OpenGL texture handle. */
	inline void init()
	{
	  if (_valid) M_throw() << "Already init()ed!";
	  glGenTextures(1, &_handle);
	  _valid = true;
	}

      public:
	/*! \brief Releases the OpenGL texture resources. */
	inline void deinit()
	{
	  if (_valid)
	    {
	      glDeleteTextures(1, &_handle);
	      _valid = false;
	    }
	}

	/*! \brief Binds the texture to the specified texture unit.
	 *
	 * \param unit The number of the texture unit to bind the
	 * texture to.
	 */
	inline void bind(int unit)
	{
	  glActiveTextureARB(GL_TEXTURE0 + unit);
	  glBindTexture(_texType, _handle);
	}

	inline void genMipmaps()
	{
	  bind(0);
	  glGenerateMipmap(_texType);
	}

	/*! \brief Sets an integer parameter of the texture.
	 *
	 * \param paramname The name of the parameter to set.
	 * \param param The value of the parameter.
	 */
	inline void parameter(GLenum paramname, GLint param)
	{ bind(0); glTexParameteri(_texType, paramname, param); }

	/*! \brief Sets a float parameter of the texture.
	 *
	 * \param paramname The name of the parameter to set.
	 * \param param The value of the parameter.
	 */
	inline void parameter(GLenum paramname, GLfloat param)
	{ bind(0); glTexParameterf(_texType, paramname, param); }

	/*! \brief Tests if the texture has been allocated. */
	inline bool isValid() const { return _valid; }

	/*! \brief Returns the OpenGL handle for the texture. */
	inline GLuint getGLHandle() 
	{ 
	  if (!_valid) M_throw() << "Texture is not initialised yet";
	  return _handle; 
	}

      protected:
      
	/*! \brief Function which returns an appropriate format
	 * parameter given a set internalformat. 
	 */
	GLenum safeFormat(GLint internalformat)
	{
	  GLenum format = GL_RGB;
	  switch (internalformat)
	    {
	    case GL_DEPTH_COMPONENT:
	    case GL_DEPTH_COMPONENT16:
	    case GL_DEPTH_COMPONENT24:
	    case GL_DEPTH_COMPONENT32:
	      format = GL_DEPTH_COMPONENT;
	    default:
	      break;
	    }

	  return format;
	}
      
	GLuint _handle;
	bool _valid;
	GLint _internalFormat;
	const GLenum _texType;
      };
    }

    /*! \brief A 1D Texture. 
     */
    class Texture1D: public detail::TextureBasic
    {
    public:
      Texture1D(): TextureBasic(GL_TEXTURE_1D) {}
      
      /*! \brief Initializes a 1D texture.
       *
       * \param width The width of the texture in pixels.
       * \param internalformat The underlying format of the texture.
       */
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
		     safeFormat(_internalFormat), GL_UNSIGNED_BYTE, NULL); 
      }

      /*! \brief Resize the texture.
       */
      inline void resize(GLint width)
      {
	//Skip identity operations
	if (width == _width) return;

	if (!_width)
	  M_throw() << "Cannot resize an uninitialised texture";

	_width = width;

	bind(0);
	glTexImage1D(_texType, 0, _internalFormat, _width,
		     0,//Border is off
		     //The following values are not used as no data is
		     //passed here, use subImage for that.
		     safeFormat(_internalFormat), GL_UNSIGNED_BYTE, NULL); 
      }

      /*! \brief Fills a section of the texture with the passed data
       *
       * \param data The pixel data to fill the texture with.
       * \param pixelformat The format of the pixel data.
       * \param xoffset The starting x position of the texture to write to.
       * \param width The amount of pixels to write.
       * \param level The texture mipmap level to write to.
       */
      inline void subImage(const std::vector<uint8_t>& data, GLenum pixelformat,
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


    /*! \brief A 2D Texture. 
     */
    class Texture2D: public detail::TextureBasic
    {
    public:
      Texture2D(): TextureBasic(GL_TEXTURE_2D) {}
      
      /*! \brief Initializes a 2D texture.
       *
       * \param width The width of the texture in pixels.
       * \param height The height of the texture in pixels.
       * \param internalformat The underlying format of the texture.
       */
      inline void init(size_t width, size_t height, GLint internalformat = GL_RGBA8)
      {
	_width = width; 
	_height = height; 
	_internalFormat = internalformat;

	TextureBasic::init();
	bind(0);
	
	glTexImage2D(_texType, 0, _internalFormat, _width, _height,
		     0, //Border is off
		     //The following values are not used as no data is
		     //passed here, use subImage for that.
		     safeFormat(_internalFormat), GL_UNSIGNED_BYTE, NULL); 
      }
      
      /*! \brief Resize the texture.
       */
      inline void resize(GLint width, GLint height)
      {
	//Skip identity operations
	if ((width == _width) && (height == _height)) return;

	if (!_width) M_throw() << "Cannot resize an uninitialised texture";

	_width = width; 
	_height = height;
	bind(0);
	glTexImage2D(_texType, 0, _internalFormat, _width, _height,
		     0, //Border is off
		     //The following values are not used as no data is
		     //passed here, use subImage for that.
		     safeFormat(_internalFormat), GL_UNSIGNED_BYTE, NULL); 	
      }

      /*! \brief Fills a section of the texture with the passed data
       *
       * \param data The pixel data to fill the texture with.
       * \param pixelformat The format of the pixel data.
       * \param xoffset The starting x position of the texture to write to.
       * \param yoffset The starting y position of the texture to write to.
       * \param width The amount of pixels to write in the x direction.
       * \param height The amount of pixels to write in the y direction.
       * \param level The texture mipmap level to write to.
       */
      inline void subImage(const std::vector<uint8_t>& data, GLenum pixelformat,
			   GLint xoffset = 0, GLint yoffset = 0, GLint width = -1, 
			   GLint height = -1, GLint level = 0)
      { 
	//If there are negative values, assume they mean to use the
	//full space
	if (width  < 0) width  = _width;
	if (height  < 0) height  = _height;
	if (xoffset < 0) M_throw() << "x offset is negative";
	if (yoffset < 0) M_throw() << "y offset is negative";
	if (xoffset + width > _width) M_throw() << "Texture write x overrun";
	if (yoffset + height > _height) M_throw() << "Texture write y overrun";

	bind(0);
	glTexSubImage2D(_texType, level, xoffset, yoffset, width, height,
			pixelformat, GL_UNSIGNED_BYTE, &data[0]);
      }

      inline void subImage(const uint8_t* data, GLenum pixelformat, GLint width, 
			   GLint height, GLint xoffset = 0, GLint yoffset = 0, GLint level = 0)
      { 
	//If there are negative values, assume they mean to use the
	//full space
	if (xoffset < 0) M_throw() << "x offset is negative";
	if (yoffset < 0) M_throw() << "y offset is negative";
	if (xoffset + width > _width) M_throw() << "Texture write x overrun";
	if (yoffset + height > _height) M_throw() << "Texture write y overrun";

	bind(0);
	glTexSubImage2D(_texType, level, xoffset, yoffset, width, height,
			pixelformat, GL_UNSIGNED_BYTE, data);
      }

      inline const GLint& getWidth() const { return _width; }
      inline const GLint& getHeight() const { return _height; }
    private:
      
      GLint _width;
      GLint _height;
    };

    /*! \brief A 3D Texture. 
     */
    class Texture3D: public detail::TextureBasic
    {
    public:
      Texture3D(): TextureBasic(GL_TEXTURE_3D) {}
      
      /*! \brief Initializes a 3D texture.
       *
       * \param width The width of the texture in pixels.
       * \param height The height of the texture in pixels.
       * \param depth The depth of the texture in pixels.
       * \param internalformat The underlying format of the texture.
       */
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
		     safeFormat(_internalFormat), GL_UNSIGNED_BYTE, NULL); 
      }

      /*! \brief Resize the texture.
       */
      inline void resize(GLint width, GLint height, GLint depth)
      {
	//Skip identity operations
	if ((width == _width) && (height == _height) && (depth == _depth))
	  return;

	if (!_width) M_throw() << "Cannot resize an uninitialised texture";

	_width = width; _height = height; _depth = depth;

	bind(0);
	glTexImage3D(_texType, 0, _internalFormat,
		     _width, _height, _depth, 
		     //Border is off
		     0,
		     //The following values are not used as no data is
		     //passed here, use subImage for that.
		     safeFormat(_internalFormat), GL_UNSIGNED_BYTE, NULL);
      }

      /*! \brief Fills a section of the texture with the passed data
       *
       * \param data The pixel data to fill the texture with.
       * \param pixelformat The format of the pixel data.
       * \param xoffset The starting x position of the texture to write to.
       * \param yoffset The starting y position of the texture to write to.
       * \param zoffset The starting z position of the texture to write to.
       * \param width The amount of pixels to write in the x direction.
       * \param height The amount of pixels to write in the y direction.
       * \param depth The amount of pixels to write in the z direction.
       * \param level The texture mipmap level to write to.
       */
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
