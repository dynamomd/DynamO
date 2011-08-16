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

#include <magnet/GL/FBO.hpp>

namespace magnet {
  namespace GL {
    /*! \brief A multisampled (anti-aliased) Frame Buffer Object.
     *
     * Multisampled FBO's use sub-pixels to render a scene at a higher
     * accuracy than is required. These sub-pixels are then averaged
     * to "smooth" the final image and remove jagged edges.
     *
     * \sa FBO
     */
    class MultisampledFBO: public FBO
    {
    public:
      /*! \brief Constructor which sets the number of subpixels per
       * pixel.
       *
       * \param samples The number of subsamples per pixel.
       */
      inline MultisampledFBO(GLsizei samples = 1):
	_samples(samples)
      {}

      inline void setSamples(GLsizei samples) { _samples = samples; }

      inline 
      virtual void init(GLsizei width, GLsizei height, GLint internalformat = GL_RGBA)
      {
	if (!GLEW_EXT_framebuffer_multisample)
	  M_throw() << "GLEW_EXT_framebuffer_multisample is not supported, cannot perform anti-aliasing";

	//Setup the RTT FBO
	FBO::init(width, height, internalformat);

	//Now create the multisampling color buffer
	glGenRenderbuffersEXT(1, &_multisampleColorBuffer);
	glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, _multisampleColorBuffer);
	glRenderbufferStorageMultisampleEXT(GL_RENDERBUFFER_EXT, _samples, internalformat, _width, _height);

	// multi sampled depth buffer
	glGenRenderbuffersEXT(1, &_multisampleDepthBuffer);
	glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, _multisampleDepthBuffer);
	glRenderbufferStorageMultisampleEXT(GL_RENDERBUFFER_EXT, _samples, GL_DEPTH_COMPONENT, _width, _height);

	//Now build the multisample FBO
	glGenFramebuffersEXT(1, &_multisampleFBO);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _multisampleFBO);
	glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, 
				     GL_RENDERBUFFER_EXT, _multisampleColorBuffer);
	glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, 
				     GL_RENDERBUFFER_EXT, _multisampleDepthBuffer);

	// check FBO status
	GLenum FBOstatus = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
	if(FBOstatus != GL_FRAMEBUFFER_COMPLETE_EXT)
	  M_throw() << "GL_FRAMEBUFFER_COMPLETE_EXT failed";
	
	// switch back to window-system-provided framebuffer
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
      }

      inline 
      virtual void init()
      {
      }
      
      inline 
      virtual void resize(GLsizei width, GLsizei height)
      {
	//If we've not been initialised, then just return
	if (!_width) return;

	//Skip identity operations
	if ((_width == width) && (_height == height)) return;

	//We don't use the resize commands, as it causes problems on
	//AMD hardware (the buffer doesn't actually resize). Instead
	//we just recreate the render buffers. This is also noted on
	//the opengl wiki.
	// http://www.opengl.org/wiki/Renderbuffer_Object
	deinit();
	init(width, height);
      }

      inline ~MultisampledFBO()
      { deinit(); }

      inline 
      virtual void attach()
      {
	if (!_width)
	  M_throw() << "Cannot attach() an uninitialised multisampledFBO";
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _multisampleFBO);
	_context->setViewport(0, 0, _width, _height);
      }

      inline 
      virtual void detach()
      {
	//First blit between the two FBO's
	glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, _multisampleFBO);
	glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER_EXT, _FBO);
	glBlitFramebufferEXT(0, 0, _width, _height, 0, 0, 
			     _width, _height, GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT, GL_NEAREST);

	glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, 0);
	glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER_EXT, 0);

	FBO::detach();
      }

      inline 
      virtual void copyto(FBO& other, GLbitfield opts = GL_COLOR_BUFFER_BIT 
			  | GL_DEPTH_BUFFER_BIT)
      {
	//First blit between the two FBO's
	glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, _multisampleFBO);
	glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER_EXT, other.getFBO());
	glBlitFramebufferEXT(0, 0, _width, _height, 0, 0, 
			     _width, _height, opts, GL_NEAREST);
	
	glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, 0);
	glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER_EXT, 0);
      }

      inline 
      virtual void deinit()
      {
	if (_width)
	  {
	    glDeleteFramebuffersEXT(1, &_multisampleFBO);
	    glDeleteRenderbuffersEXT(1, &_multisampleDepthBuffer);
	    glDeleteRenderbuffersEXT(1, &_multisampleColorBuffer);
	  }

	FBO::deinit();
      }
      
      /*! \brief Returns the number of sub-samples supported by the OpenGL implementation.
       *
       * This function returns 1, if multisampling is not supported.
       */
      static GLint getSupportedSamples()
      {
	if (!GLEW_EXT_framebuffer_multisample) return 1;

	GLint maxSamples;
	glGetIntegerv(GL_MAX_SAMPLES, &maxSamples);
	return maxSamples;
      }

    private:
      GLuint _multisampleFBO;
      GLuint _multisampleColorBuffer;
      GLuint _multisampleDepthBuffer;
      GLsizei _samples;
    };
  }
}
