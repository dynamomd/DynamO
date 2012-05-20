/*    dynamo:- Event driven molecular dynamics simulator 
 *    http://www.dynamomd.org
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

      /*! \brief Set the number of samples to be used by the multisampling buffers
       */
      inline void setSamples(GLsizei samples) 
      {
	_samples = samples;
	_validated = false;
      }

      inline virtual void init()
      {
	if (!GLEW_EXT_framebuffer_multisample)
	  M_throw() << "GLEW_EXT_framebuffer_multisample is not supported, cannot perform anti-aliasing";

	FBO::init();
	glGenFramebuffersEXT(1, &_multisampleFBO);
	_colorRenderBuffers.resize(_colorTextures.size());
      }

      inline ~MultisampledFBO() { deinit(); }

      inline 
      virtual void attach()
      {
	if (!_context)
	  M_throw() << "Cannot attach() an uninitialised multisampledFBO";
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _multisampleFBO);
	_context->setViewport(0, 0, getWidth(), getHeight());

	std::vector<GLenum> states(_colorTextures.size(), GL_NONE);
	for (size_t attachment(0); attachment < _colorTextures.size(); ++attachment)
	  if (_colorTextures[attachment])
	    states[attachment] = GL_COLOR_ATTACHMENT0_EXT + attachment;
	glDrawBuffers(states.size(), &states[0]);
      }

      inline 
      virtual void detach()
      {
	validate();
	if (!_context)
	  M_throw() << "Cannot detach() an uninitialised multisampledFBO";
	//First blit between the two FBO's
	glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, _multisampleFBO);
	glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER_EXT, _FBO);
	glBlitFramebufferEXT(0, 0, getWidth(), getHeight(), 0, 0, 
			     getWidth(), getHeight(), 
			     GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT, GL_NEAREST);

	glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, 0);
	glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER_EXT, 0);

	FBO::detach();
      }

      inline 
      virtual void copyto(FBO& other, GLbitfield opts = GL_COLOR_BUFFER_BIT 
			  | GL_DEPTH_BUFFER_BIT)
      {
	validate();
	//First blit between the two FBO's
	glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, _multisampleFBO);
	glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER_EXT, other.getFBO());
	glBlitFramebufferEXT(0, 0, getWidth(), getHeight(), 0, 0, 
			     getWidth(), getHeight(), opts, GL_NEAREST);
	
	glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, 0);
	glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER_EXT, 0);
      }

      virtual void deinit()
      {
	_colorRenderBuffers.clear();
	_depthRenderBuffer.deinit();
	if (_context)
	  glDeleteFramebuffersEXT(1, &_multisampleFBO);

	FBO::deinit();
      }
      
      /*! \brief Returns the number of sub-samples supported by the OpenGL implementation.
       *
       * This function returns 1, if multisampling is not supported.
       */
      static GLint getSupportedSamples()
      {
	if (!GLEW_EXT_framebuffer_multisample) return 1;
	return detail::glGet<GL_MAX_SAMPLES>();
      }

    private:
      /*! \brief Class to represent and manage a single (possibly
          multisampling) render buffer.
       */
      class RenderBuffer
      {
      public:
	/*! \brief Constructor
	 */
	RenderBuffer(): _valid(false) {}

	/*! \brief Destructor
	 */
	~RenderBuffer() { deinit(); }

	/*! \brief Initialises the OpenGL resources for this render buffer.
	  
	  \param width The width of the render buffer in pixels.
	  \param height The height of the render buffer in pixels.
	  \param internalformat The pixel format of the buffer (e.g., GL_RGBA).
	  \param samples The number of pixel sub-samples (0 disables multisampling).
	 */
	virtual void init(GLsizei width, GLsizei height, GLint internalformat, GLsizei samples)
	{
	  deinit();
	  glGenRenderbuffersEXT(1, &_buf);
	  bind();
	  glRenderbufferStorageMultisampleEXT(GL_RENDERBUFFER_EXT, samples, internalformat, width, height);
	  _valid = true;
	}

	/*! \brief Release the resources obtained by this buffer.
	 */
	virtual void deinit()
	{
	  if (_valid) glDeleteRenderbuffersEXT(1, &_buf);
	  _valid = false;
	}

	/*! \brief Bind the renderbuffer to the OpenGL state.
	 */
	void bind() { glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, _buf); }

	void attach(GLenum attachment)
	{ glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, attachment, GL_RENDERBUFFER_EXT, _buf); }

      protected:
	bool _valid;
	GLuint _buf;
      };

      virtual void validate()
      {
	if (!_context)
	  M_throw() << "Cannot attach() an uninitialised FBO";

	//We let the underlying FBO validate first as it will verify
	//the bound texture formats etc.
	bool alreadyvalidated = _validated;
	FBO::validate();

	//Now rebuild the multisampled attachments
	if(!alreadyvalidated)
	  {
	    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _multisampleFBO);

	    if (_depthTexture)
	      {
		_depthRenderBuffer.init(_depthTexture->getWidth(), _depthTexture->getHeight(),
					_depthTexture->getInternalFormat(), _samples);
		_depthRenderBuffer.attach(GL_DEPTH_ATTACHMENT_EXT);
	      }
	    else
	      {
		glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, 0);
		_depthRenderBuffer.deinit();
	      }
		
	    for (size_t attachment(0); attachment < _colorTextures.size(); ++attachment)
	      if (_colorTextures[attachment])
		{
		  _colorRenderBuffers[attachment].init(_colorTextures[attachment]->getWidth(), 
						       _colorTextures[attachment]->getHeight(),
						       _colorTextures[attachment]->getInternalFormat(), _samples);
		  _colorRenderBuffers[attachment].attach(GL_COLOR_ATTACHMENT0_EXT + attachment);
		}
	      else
		{
		  glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT + attachment,
					       GL_RENDERBUFFER_EXT, 0);
		  _colorRenderBuffers[attachment].deinit();
		}

	    // check FBO status
	    GLenum FBOstatus = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
	    if(FBOstatus != GL_FRAMEBUFFER_COMPLETE_EXT)
	      M_throw() << "GL_FRAMEBUFFER_COMPLETE_EXT failed";
	  }
      }

      GLuint _multisampleFBO;
      std::vector<RenderBuffer> _colorRenderBuffers;
      RenderBuffer _depthRenderBuffer;
      GLsizei _samples;
    };
  }
}
