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

namespace magnet {
  namespace GL {    
    class FBO
    {
    public:
      inline FBO():
	_width(0),
	_height(0)
      {}

      inline void init(GLsizei width, GLsizei height, GLint internalformat = GL_RGBA, 
		       GLenum format = GL_RGBA, GLenum type = GL_UNSIGNED_BYTE)
      {
	if (!GLEW_EXT_framebuffer_object)
	  M_throw() << "GL_EXT_framebuffer_object extension is not supported! Cannot do offscreen rendering!";

	if (_width || _height) 
	  M_throw() << "FBO has already been initialised!";

	if (!width || !height)
	  M_throw() << "Cannot initialise an FBO with a width or height == 0!";

	_internalformat = internalformat;
	_format = format;	
	_type = type;		
	_width = width;
	_height = height;

	//Build depth buffer
	glGenRenderbuffersEXT(1, &_depthBuffer);
	glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, _depthBuffer);
	glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT, _width, _height);

	//Build color texture
	glGenTextures(1, &_colorTexture);
	glBindTexture(GL_TEXTURE_2D, _colorTexture);
	glTexImage2D(GL_TEXTURE_2D, 0, _internalformat, _width, _height, 0, _format, _type, NULL);

	//Build the FBO
	glGenFramebuffersEXT(1, &_FBO);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _FBO);

	//Bind the depth buffer
	glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, 
				     _depthBuffer);

	//Bind the texture
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, _colorTexture, 0);

	// check FBO status
	GLenum FBOstatus = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
	if(FBOstatus != GL_FRAMEBUFFER_COMPLETE_EXT)
	  M_throw() << "GL_FRAMEBUFFER_COMPLETE_EXT failed";
	
	// switch back to window-system-provided framebuffer
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);

      }
      
      inline void resize(GLsizei width, GLsizei height)
      {
	//If we've not been initialised, then just return
	if (!_width) return;

	_width = width;
	_height = height;

	glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, _depthBuffer);
	glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT, _width, _height);
	
	glBindTexture(GL_TEXTURE_2D, _colorTexture);
	glTexImage2D(GL_TEXTURE_2D, 0, _internalformat, _width, _height, 0, _format, _type, NULL);
      }

      inline void blitToScreen(GLsizei screenwidth, GLsizei screenheight)
      {
	if (!GLEW_EXT_framebuffer_blit)
	  M_throw() << "The GL_EXT_framebuffer_blit extension is not supported! Cannot blit!";
	glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, _FBO);
	glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER_EXT, 0);
	glBlitFramebufferEXT(0, 0, _width, _height, 0, 0, screenwidth, screenheight, 
			     GL_COLOR_BUFFER_BIT, GL_NEAREST);
	glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, 0);
      }

      inline ~FBO()
      {
	if (_width)
	  {
	    glDeleteFramebuffersEXT(1, &_FBO);
	    glDeleteTextures(1, &_colorTexture);
	    glDeleteRenderbuffersEXT(1, &_depthBuffer);
	  }
      }

      inline void attach()
      {
	glPushAttrib(GL_VIEWPORT_BIT);
	glViewport(0, 0, _width, _height);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _FBO);
      }

      inline void detach()
      {
	glPopAttrib();
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
      }
	

      inline GLuint getFBO() { return _FBO; }
      inline GLuint getColorTexture() { return _colorTexture; }
      inline GLsizei getWidth() { return _width; }
      inline GLsizei getHeight() { return _height; }

    protected:
      GLuint _FBO, _colorTexture, _depthBuffer;
      GLsizei _width, _height;
      GLint _internalformat;
      GLenum _format;	
      GLenum _type;		

    };
  }
}
