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

namespace magnet {
  namespace GL {    
    class FBO
    {
    public:
      inline FBO():
	_width(0),
	_height(0)
      {}

      inline 
      virtual void init(GLsizei width, GLsizei height, GLint internalformat = GL_RGBA)
      {
	if (!GLEW_EXT_framebuffer_object)
	  M_throw() << "GL_EXT_framebuffer_object extension is not supported! Cannot do offscreen rendering!";

	if (_width || _height) 
	  M_throw() << "FBO has already been initialised!";

	if (!width || !height)
	  M_throw() << "Cannot initialise an FBO with a width or height == 0!";

	_internalformat = internalformat;
	_width = width;
	_height = height;

	//Build depth buffer
	glGenTextures(1, &_depthTexture);
	glBindTexture(GL_TEXTURE_2D, _depthTexture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, _width, _height, 0,
		     GL_DEPTH_COMPONENT, GL_UNSIGNED_INT, NULL);

	//Where to put the information
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_NONE);
	//////////////////////////////////////////////////////////////////////////

	//Build color texture
	glGenTextures(1, &_colorTexture);
	glBindTexture(GL_TEXTURE_2D, _colorTexture);
	glTexImage2D(GL_TEXTURE_2D, 0, _internalformat, _width, _height, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	glGenFramebuffersEXT(1, &_FBO);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _FBO);

	//Bind the depth texture
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, 
				  GL_TEXTURE_2D, _depthTexture, 0);

	//Bind the color texture
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, _colorTexture, 0);

	// check FBO status
	GLenum FBOstatus = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
	if(FBOstatus != GL_FRAMEBUFFER_COMPLETE_EXT)
	  M_throw() << "GL_FRAMEBUFFER_COMPLETE_EXT failed";
	
	// switch back to window-system-provided framebuffer
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
      }
      
      inline 
      virtual void resize(GLsizei width, GLsizei height)
      {
	//If we've not been initialised, then just return
	if (!_width) return;
	
	//Skip identity operations
	if ((_width == width) && (_height == height)) return;

	_width = width;
	_height = height;

	glBindTexture(GL_TEXTURE_2D, _depthTexture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, _width, _height, 0,
		     GL_DEPTH_COMPONENT, GL_UNSIGNED_INT, NULL);

	glBindTexture(GL_TEXTURE_2D, _colorTexture);
	glTexImage2D(GL_TEXTURE_2D, 0, _internalformat, _width, _height, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
      }

      inline void blitToScreen(GLsizei screenwidth, GLsizei screenheight)
      {
	if (!GLEW_EXT_framebuffer_blit)
	  M_throw() << "The GL_EXT_framebuffer_blit extension is not supported! Cannot blit!";
	glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, _FBO);
	glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER_EXT, 0);
	glBlitFramebufferEXT(0, 0, _width, _height, 0, 0, screenwidth, screenheight, 
			     GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT, GL_NEAREST);
	glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, 0);
      }

      inline 
      virtual ~FBO()
      { deinit(); }

      inline 
      virtual void deinit()
      {
	if (_width)
	  {
	    glDeleteFramebuffersEXT(1, &_FBO);
	    glDeleteTextures(1, &_colorTexture);
	    glDeleteTextures(1, &_depthTexture);
	  }

	_width = 0;
	_height = 0;
      }

      inline 
      virtual void attach()
      {
	glPushAttrib(GL_VIEWPORT_BIT);
	glViewport(0, 0, _width, _height);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _FBO);
      }

      inline 
      virtual void detach()
      {
	glPopAttrib();
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
      }

      inline 
      virtual void copyto(FBO& other, GLbitfield opts = GL_COLOR_BUFFER_BIT 
			  | GL_DEPTH_BUFFER_BIT)
      {
	//First blit between the two FBO's
	glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, _FBO);
	glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER_EXT, other._FBO);
	glBlitFramebufferEXT(0, 0, _width, _height, 0, 0, 
			     _width, _height, opts, GL_NEAREST);
	
	glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, 0);
	glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER_EXT, 0);
      }

      
      inline GLuint getFBO() { return _FBO; }
      inline GLuint getColorTexture() { return _colorTexture; }
      inline GLuint getDepthTexture() { return _depthTexture; }
      inline GLsizei getWidth() { return _width; }
      inline GLsizei getHeight() { return _height; }

    protected:
      GLuint _FBO, _colorTexture, _depthTexture;
      GLsizei _width, _height;
      GLint _internalformat;
    };
  }
}
