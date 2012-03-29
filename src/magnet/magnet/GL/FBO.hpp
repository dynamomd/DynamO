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
     
      Frame buffer objects are "virtual screens" which can be drawn
      to, but the output is captured to bound textures instead of the
      users screen.
     
      This framebuffer wrapper uses a validate-on-attachment
      methodology like the underlying OpenGL FBO. This means that you
      initialise the FBO, attach buffers to its attachment points and
      when you call attach() it validates the configuration.

      \sa MultisampledFBO
     */
    class FBO
    {
    public:
      /*! \brief Default constructor.
       * 
       * The FBO is unusable at this point and must be first \ref
       * init() ialized.
       */
      inline FBO(): _validated(false) {}

      /*! \brief Initializes the FBO.
	
	This function does not attach any textures to the FBO, you
	must attach at least one texture using \ref attachColorTexture() or \ref \attachDepthTexture() 
	before attempting to \ref attach() this FBO.
       */
      inline virtual void init()
      {
	if (_context)
	  M_throw() << "FBO has already been initialised!";
	
	_context = Context::getContext();

	glGenFramebuffersEXT(1, &_FBO);
	detail::errorCheck();
	//Here we only allocate enough texture pointers for the maximum
	//allowed drawable buffers!
	_colorTextures.resize(detail::glGet<GL_MAX_DRAW_BUFFERS>());
	_validated = false;
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
      inline void resize(const GLsizei width, const GLsizei height)
      {
	if (!_context)
	  M_throw() << "Cannot resize an uninitialized FBO";

	//We don't check if the FBO is validated as textures might be
	//shared between multiple FBO's. We should resize if *any* of
	//the bound textures of this FBO do not match the new size.
	{
	  //Check if the resize is a no_op
	  size_t no_op = true;
	  for (size_t attachment(0); attachment < _colorTextures.size(); ++attachment)
	    if (_colorTextures[attachment])
	      if ((_colorTextures[attachment]->getWidth() != width)
		  || (_colorTextures[attachment]->getHeight() != height))
		{ no_op = false; break; }

	  if (_depthTexture)
	    if ((_depthTexture->getWidth() != width) 
		|| (_depthTexture->getHeight() != height))
	      no_op = false;
	  
	  if (no_op) return;
	}

	std::vector<std::tr1::shared_ptr<Texture2D> > colorTextures = _colorTextures;
	std::tr1::shared_ptr<Texture2D> depthTexture = _depthTexture;

	deinit();
	init();

	for (size_t attachment(0); attachment < colorTextures.size(); ++attachment)
	  if (colorTextures[attachment]) 
	    { 
	      colorTextures[attachment]->resize(width, height); 
	      attachTexture(colorTextures[attachment], attachment);
	    }

	if (depthTexture)
	  {
	    depthTexture->resize(width, height);
	    attachTexture(depthTexture);
	  }
      }

      /*! \brief Renders the contents of the FBO to the real screen FBO.
       * 
       * \param screenwidth The width of the screen in pixels.
       * \param screenheight The height of the screen in pixels.
       */
      inline void blitToScreen(GLsizei screenwidth, GLsizei screenheight, 
			       GLsizei screenx = 0, GLsizei screeny = 0,
			       GLenum filter = GL_NEAREST,
			       GLbitfield mask = GL_COLOR_BUFFER_BIT)
      {
	validate();

	if (!GLEW_EXT_framebuffer_blit)
	  M_throw() << "The GL_EXT_framebuffer_blit extension is not supported!"
		    << " Cannot blit!";
	glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, _FBO);
	detail::errorCheck();
	glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER_EXT, 0);
	detail::errorCheck();
	glBlitFramebufferEXT(0, 0, getWidth(), getHeight(), screenx, screeny, 
			     screenwidth + screenx, screenheight + screeny, 
			     mask, filter);
	detail::errorCheck();
	glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, 0);
	detail::errorCheck();
      }

      inline virtual ~FBO() { deinit(); }

      /*! \brief Releases the OpenGL resources of this FBO. */
      inline 
      virtual void deinit()
      {
	_colorTextures.clear();
	_depthTexture.reset();
	
	if (_context)
	  {
	    glDeleteFramebuffersEXT(1, &_FBO);
	    detail::errorCheck();
	  }

	_context.reset();
	_validated = false;
      }

      /*! \brief Attaches this FBO as the current render target. */
      inline virtual void attach()
      {
	validate();
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _FBO);
	detail::errorCheck();
	_context->setViewport(0, 0, getWidth(), getHeight());

	std::vector<GLenum> states(_colorTextures.size(), GL_NONE);
	for (size_t attachment(0); attachment < _colorTextures.size(); ++attachment)
	  if (_colorTextures[attachment])
	    states[attachment] = GL_COLOR_ATTACHMENT0_EXT + attachment;

	glDrawBuffers(states.size(), &states[0]);
	detail::errorCheck();
      }

      /*! \brief Restores the screen FBO as the current render target. */
      inline 
      virtual void detach()
      {
	if (!_context)
	  M_throw() << "Cannot detach() an uninitialised FBO";

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
	detail::errorCheck();
	//No need to reset the draw states, glDrawBuffers is not valid
	//for the default framebuffer.
      }

      inline void attachTexture(std::tr1::shared_ptr<Texture2D> tex, size_t i = 0)
      {
	if (!_context)
	  M_throw() << "Cannot attach textures to an uninitialised FBO";
	
	switch (tex->getInternalFormat())
	  {
	  case GL_DEPTH24_STENCIL8:
	  case GL_DEPTH32F_STENCIL8:
	  case GL_DEPTH_COMPONENT24:
	  case GL_DEPTH_COMPONENT32:
	  case GL_DEPTH_COMPONENT32F:
	  case GL_DEPTH_COMPONENT:
	    if (i)
	      M_throw() << "Texture attachment point out of range";
	    _depthTexture = tex;
	    break;
	  default:
	    if (i >= _colorTextures.size())
	      M_throw() << "Texture attachment point Out of range";
	    _colorTextures[i] = tex;
	    break;
	  }
	
	_validated = false;
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
	validate();
	other.validate();
	//First blit between the two FBO's
	glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, _FBO);
	detail::errorCheck();
	glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER_EXT, other._FBO);
	detail::errorCheck();

	glBlitFramebufferEXT(0, 0,       getWidth(),       getHeight(), 
			     0, 0, other.getWidth(), other.getHeight(),
			     opts, GL_NEAREST);
	detail::errorCheck();

	glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, 0);
	detail::errorCheck();
	glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER_EXT, 0);
	detail::errorCheck();
      }

      
      /*! \brief Fetch the underlying OpenGL handle. */
      inline GLuint getFBO() { return _FBO; }

      /*! \brief Fetch the texture bound to the color buffer. */
      inline std::tr1::shared_ptr<Texture2D>& getColorTexture(const size_t ID = 0)
      { 
	if (ID >= _colorTextures.size())
	  M_throw() << "Out of range";

	if (!_colorTextures[ID]) 
	  M_throw() << "Cannot fetch the color texture " << ID << " as the FBO has none bound";

	return _colorTextures[ID];
      }

      /*! \brief Fetch the texture bound to the depth buffer. */
      inline std::tr1::shared_ptr<Texture2D>& getDepthTexture()
      { 
	if (!_depthTexture) 
	  M_throw() << "Cannot fetch the depth texture as the FBO has none bound";
	return _depthTexture; 
      }

      /*! \brief Fetch the width of the FBO in pixels. */
      inline GLsizei getWidth() 
      {
	validate(); //Check the format of the FBO is consistent!
	//Find the first bound texture and return its dimensions
	if (_depthTexture) return _depthTexture->getWidth();

	for (std::vector<std::tr1::shared_ptr<Texture2D> >::iterator iPtr = _colorTextures.begin();
	     iPtr != _colorTextures.end(); ++iPtr)
	  if (*iPtr) return (*iPtr)->getWidth();
	
	M_throw() << "Cannot query the width of a FBO without any bound textures";
      }

      /*! \brief Fetch the height of the FBO in pixels. */
      inline GLsizei getHeight() 
      {
	validate(); //Check the format of the FBO is consistent!
	//Find the first bound texture and return its dimensions
	if (_depthTexture) return _depthTexture->getHeight();

	for (std::vector<std::tr1::shared_ptr<Texture2D> >::iterator iPtr = _colorTextures.begin();
	     iPtr != _colorTextures.end(); ++iPtr)
	  if (*iPtr) return (*iPtr)->getHeight();
	
	M_throw() << "Cannot query the height of a FBO without any bound textures";
      }

      Context& getContext()
      { 
	if (!_context)
	  M_throw() << "Cannot get an FBO's context if it is uninitialized";
	return *_context;
      }

      void checkStatus()
      {
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _FBO);
	detail::errorCheck();

	// check FBO status
	GLenum FBOstatus = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
	switch (FBOstatus)
	  {
	  case GL_FRAMEBUFFER_UNDEFINED: 
	    M_throw() << "Failed to create FrameBufferObject: GL_FRAMEBUFFER_UNDEFINED";
	  case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT: 
	    M_throw() << "Failed to create FrameBufferObject: GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT";
	  case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT: 
	    M_throw() << "Failed to create FrameBufferObject: GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT";
	  case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER: 
	    M_throw() << "Failed to create FrameBufferObject: GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER";
	  case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER: 
	    M_throw() << "Failed to create FrameBufferObject: GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER";
	  case GL_FRAMEBUFFER_UNSUPPORTED: 
	    M_throw() << "Failed to create FrameBufferObject: GL_FRAMEBUFFER_UNSUPPORTED";
	  case GL_FRAMEBUFFER_INCOMPLETE_MULTISAMPLE: 
	    M_throw() << "Failed to create FrameBufferObject: GL_FRAMEBUFFER_INCOMPLETE_MULTISAMPLE";
	  case GL_FRAMEBUFFER_INCOMPLETE_LAYER_TARGETS: 
	    M_throw() << "Failed to create FrameBufferObject: GL_FRAMEBUFFER_INCOMPLETE_LAYER_TARGETS";
	  default: 
	    M_throw() << "Failed to create FrameBufferObject: Unknown error code = " << FBOstatus << " ";
	  case GL_FRAMEBUFFER_COMPLETE_EXT: 
	    break;
	  }
      }


    protected:
      /*! \brief Validate the current state of the FBO and raise an exception if there is an error.
       */
      virtual void validate()
      {
	if (!_context)
	  M_throw() << "Cannot attach() an uninitialised FBO";

	if(!_validated)
	  {
	    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _FBO);
	    detail::errorCheck();

	    //Bind the textures, or unbind the unbound textures ready for the completeness test
	    if (_depthTexture)
	      {
		if ((_depthTexture->getInternalFormat() == GL_DEPTH24_STENCIL8)
		    || (_depthTexture->getInternalFormat() == GL_DEPTH32F_STENCIL8))
		  {
		    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_STENCIL_ATTACHMENT, 
					      _depthTexture->getGLType(), _depthTexture->getGLHandle(), 0);
		    detail::errorCheck();
		  }
		else
		  {
		    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, 
					      _depthTexture->getGLType(), _depthTexture->getGLHandle(), 0);
		    detail::errorCheck();
		  }
	      }
	    else
	      {
		glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, 
					  GL_TEXTURE_2D, 0, 0);
		detail::errorCheck();
	      }

	    
	    std::vector<GLenum> states(_colorTextures.size(), GL_NONE);

	    for (size_t attachment(0); attachment < _colorTextures.size(); ++attachment)
	      if (_colorTextures[attachment])
		{
		  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT + attachment, 
					    _colorTextures[attachment]->getGLType(), 
					    _colorTextures[attachment]->getGLHandle(), 0);
		  detail::errorCheck();
		  states[attachment] = GL_COLOR_ATTACHMENT0_EXT + attachment;
		}
	      else
		{
		  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT + attachment, GL_TEXTURE_2D, 
					    0, 0);
		  detail::errorCheck();
		}

	    glDrawBuffers(states.size(), &states[0]);
	    detail::errorCheck();
	    checkStatus();
	    _validated = true;
	  }
      }      

      Context::ContextPtr _context;
      std::vector<std::tr1::shared_ptr<Texture2D> > _colorTextures;
      std::tr1::shared_ptr<Texture2D> _depthTexture;
      GLuint _FBO;
      bool _validated;
    };
  }
}
