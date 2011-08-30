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
#include <magnet/GL/texture.hpp>
#include <tr1/memory>

namespace magnet {
  namespace GL {   
    /*! \brief A Frame Buffer Object.
     *
     * Frame buffer objects are "virtual screens" which can be drawn
     * to, but the output is captured by bound textures instead of the
     * real user screen.
     *
     * \sa MultisampledFBO
     */
    class FBO
    {
    public:
      /*! \brief Default constructor.
       * 
       * The FBO is unusable at this point and must be first \ref
       * init() ialized.
       */
      inline FBO():_context(NULL), _width(0), _height(0) {}

      /*! \brief Initializes the FBO
       * 
       * \param width The width of the FBO in pixels.
       * \param height The height of the FBO in pixels.
       * \param internalformat The pixel type stored by the FBO.
       */
      inline 
      virtual void init(GLsizei width, GLsizei height, GLint internalformat = GL_RGBA)
      {
	deinit();
	if (!GLEW_EXT_framebuffer_object)
	  M_throw() << "GL_EXT_framebuffer_object extension is not supported! Cannot do offscreen rendering!";

	//Build depth buffer
	std::tr1::shared_ptr<Texture2D> depthTexture(new Texture2D);
	//We don't force GL_DEPTH_COMPONENT24 as it is likely you get
	//the best precision anyway
	depthTexture->init(width, height, GL_DEPTH_COMPONENT);
	//You must select GL_NEAREST for depth data, as GL_LINEAR
	//converts the value to 8bit for interpolation (on NVidia).
	depthTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	depthTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	depthTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	depthTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	depthTexture->parameter(GL_TEXTURE_COMPARE_MODE, GL_NONE);
	//////////////////////////////////////////////////////////////////////////

	//Build color texture
	std::tr1::shared_ptr<Texture2D> colorTexture(new Texture2D);
	colorTexture->init(width, height, internalformat);
	colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	colorTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	colorTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	colorTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	init(colorTexture, depthTexture);
      }

      /*! \brief Initializes the FBO using the passed color and depth buffers
       * 
       * \param width The width of the FBO in pixels.
       * \param height The height of the FBO in pixels.
       * \param internalformat The pixel type stored by the FBO.
       */
      inline 
      virtual void init(std::tr1::shared_ptr<Texture2D> colorTexture,
			std::tr1::shared_ptr<Texture2D> depthTexture)
      {
	if (_width || _height)
	  M_throw() << "FBO has already been initialised!";

	if (!colorTexture && !depthTexture)
	  M_throw() << "You must provide at least a color or a depth texture to the FBO";

	if (colorTexture && depthTexture)
	  if ((colorTexture->getWidth() != depthTexture->getWidth())
	      || (colorTexture->getHeight() != depthTexture->getHeight()))
	    M_throw() << "color and depth texture size mismatch";
	
	_context = &Context::getContext();

	glGenFramebuffersEXT(1, &_FBO);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _FBO);

	_colorTextures.resize(detail::glGet<GL_MAX_COLOR_ATTACHMENTS_EXT>());
	_colorTextures[0] = colorTexture;
	_depthTexture = depthTexture;

	if (_colorTextures[0])
	  {
	    _width = _colorTextures[0]->getWidth();
	    _height = _colorTextures[0]->getHeight();
	  }
	else
	  {
	    _width = _depthTexture->getWidth();
	    _height = _depthTexture->getHeight();
	  }
	  		
	//Bind the depth texture
	if (_depthTexture)
	  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, 
				    GL_TEXTURE_2D, _depthTexture->getGLHandle(), 0);
	
	//Bind the color texture
	if (_colorTextures[0])
	  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, 
				    _colorTextures[0]->getGLHandle(), 0);
	
	// check FBO status
	GLenum FBOstatus = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
	switch (FBOstatus)
	  {
	  case GL_FRAMEBUFFER_UNDEFINED: M_throw() << "Failed to create FrameBufferObject: GL_FRAMEBUFFER_UNDEFINED";
	  case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT: M_throw() << "Failed to create FrameBufferObject: GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT";
	  case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT: M_throw() << "Failed to create FrameBufferObject: GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT";
	  case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER: M_throw() << "Failed to create FrameBufferObject: GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER";
	  case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER: M_throw() << "Failed to create FrameBufferObject: GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER";
	  case GL_FRAMEBUFFER_UNSUPPORTED: M_throw() << "Failed to create FrameBufferObject: GL_FRAMEBUFFER_UNSUPPORTED";
	  case GL_FRAMEBUFFER_INCOMPLETE_MULTISAMPLE: M_throw() << "Failed to create FrameBufferObject: GL_FRAMEBUFFER_INCOMPLETE_MULTISAMPLE";
	  case GL_FRAMEBUFFER_INCOMPLETE_LAYER_TARGETS: M_throw() << "Failed to create FrameBufferObject: GL_FRAMEBUFFER_INCOMPLETE_LAYER_TARGETS";
	  default: M_throw() << "Failed to create FrameBufferObject: Unkown error code";
	  case GL_FRAMEBUFFER_COMPLETE_EXT: break;
	  }
      }
      
      /*! \brief Resizes the FBO
       * 
       * If the FBO has been initialized, this will resize the FBO by
       * first deinit()ializing the FBO and re-init()ializing with the
       * correct size. Just resizing the textures does not play well
       * with all GPU's.
       *
       * \param width The new width of the FBO in pixels.
       * \param height The new height of the FBO in pixels.
       */
      inline 
      virtual void resize(GLsizei width, GLsizei height)
      {
	if (!_width)
	  M_throw() << "Cannot resize an uninitialized FBO";
	
	//Skip identity operations
	if ((_width == width) && (_height == height)) return;
	
	std::vector<std::tr1::shared_ptr<Texture2D> > colorTextures = _colorTextures;
	std::tr1::shared_ptr<Texture2D> depthTexture = _depthTexture;

	deinit();
	for (std::vector<std::tr1::shared_ptr<Texture2D> >::iterator iPtr = _colorTextures.begin();
	     iPtr != _colorTextures.end(); ++iPtr)
	  if (*iPtr)
	    (*iPtr)->resize(width, height);

	if (depthTexture)
	  depthTexture->resize(width, height);

	init(colorTextures[0], depthTexture);
      }

      /*! \brief Renders the contents of the FBO to the real screen FBO.
       * 
       * \param screenwidth The width of the screen in pixels.
       * \param screenheight The height of the screen in pixels.
       */
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

      inline virtual ~FBO() { deinit(); }

      /*! \brief Releases the OpenGL resources of this FBO. */
      inline 
      virtual void deinit()
      {
	_colorTextures.clear();
	_depthTexture.reset();
	
	if (_width)
	  glDeleteFramebuffersEXT(1, &_FBO);

	_width = 0;
	_height = 0;
	_context = NULL;
      }

      /*! \brief Attaches this FBO as the current render target. */
      inline 
      virtual void attach()
      {
	if (!_width)
	  M_throw() << "Cannot attach() an uninitialised FBO";
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _FBO);
	_context->setViewport(0, 0, _width, _height);
      }

      /*! \brief Restores the screen FBO as the current render target. */
      inline 
      virtual void detach()
      {
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
      }

      /*! \brief Copies the contents of this FBO to another.
       *
       * \param other The FBO to copy this FBO's contents to.
       * \param opts Bit flags marking the channels of the FBO to copy.
       */
      inline 
      virtual void copyto(FBO& other, GLbitfield opts = GL_COLOR_BUFFER_BIT 
			  | GL_DEPTH_BUFFER_BIT)
      {
	//First blit between the two FBO's
	glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, _FBO);
	glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER_EXT, other._FBO);
	glBlitFramebufferEXT(0, 0, _width, _height, 0, 0, 
			     other._width, other._height, opts, GL_NEAREST);
	
	glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, 0);
	glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER_EXT, 0);
      }

      
      /*! \brief Fetch the underlying OpenGL handle. */
      inline GLuint getFBO() { return _FBO; }

      /*! \brief Fetch the texture bound to the color buffer. */
      inline Texture2D& getColorTexture(const size_t ID = 0) 
      { 
	if (ID >= _colorTextures.size())
	  M_throw() << "Out of range";

	if (!_colorTextures[ID]) 
	  M_throw() << "Cannot fetch the color texture " << ID << " if the FBO has none bound";
	return *_colorTextures[ID]; 
      }

      /*! \brief Fetch the texture bound to the depth buffer. */
      inline Texture2D& getDepthTexture() 
      { 
	if (!_depthTexture) M_throw() << "Cannot fetch the color texture if the FBO has none bound";
	return *_depthTexture; 
      }

      /*! \brief Fetch the width of the FBO in pixels. */
      inline GLsizei getWidth() { return _width; }

      /*! \brief Fetch the height of the FBO in pixels. */
      inline GLsizei getHeight() { return _height; }

      Context& getContext()
      { 
	if (!_width) 
	  M_throw() << "Cannot get an FBO's context if it is uninitialized";
	return *_context;
      }

    protected:
      Context* _context;
      std::vector<std::tr1::shared_ptr<Texture2D> > _colorTextures;
      std::tr1::shared_ptr<Texture2D> _depthTexture;
      GLuint _FBO;
      GLsizei _width, _height;
    };
  }
}
