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
      inline FBO():_context(NULL), _validated(false) {}

      /*! \brief Initializes the FBO.
	
	This function does not attach any textures to the FBO, you
	must attach at least one texture using \ref attachColorTexture() or \ref \attachDepthTexture() 
	before attempting to \ref attach() this FBO.
       */
      inline virtual void init()
      {
	if (_context)
	  M_throw() << "FBO has already been initialised!";
	
	_context = &Context::getContext();

	glGenFramebuffersEXT(1, &_FBO);
	//Here we only allocate enough textures for the maximum
	//allowed drawable buffers!
	//Any more would just be silly and confuse things!
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
      inline void resize(GLsizei width, GLsizei height)
      {
	if (!_context)
	  M_throw() << "Cannot resize an uninitialized FBO";
	
	//Skip identity operations
	if ((width == getWidth()) && (height == getHeight())) return;
	
	std::vector<std::tr1::shared_ptr<Texture2D> > colorTextures = _colorTextures;
	std::tr1::shared_ptr<Texture2D> depthTexture = _depthTexture;

	deinit();
	init();

	for (size_t attachment(0); attachment < colorTextures.size(); ++attachment)
	  if (colorTextures[attachment]) 
	    { 
	      colorTextures[attachment]->resize(width, height); 
	      attachColorTexture(colorTextures[attachment], attachment);
	    }

	if (depthTexture)
	  {
	    depthTexture->resize(width, height);
	    attachDepthTexture(depthTexture);
	  }
      }

      /*! \brief Renders the contents of the FBO to the real screen FBO.
       * 
       * \param screenwidth The width of the screen in pixels.
       * \param screenheight The height of the screen in pixels.
       */
      inline void blitToScreen(GLsizei screenwidth, GLsizei screenheight)
      {
	validate();

	if (!GLEW_EXT_framebuffer_blit)
	  M_throw() << "The GL_EXT_framebuffer_blit extension is not supported! Cannot blit!";
	glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, _FBO);
	glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER_EXT, 0);
	glBlitFramebufferEXT(0, 0, getWidth(), getHeight(), 0, 0, screenwidth, screenheight, 
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
	
	if (_context)
	  glDeleteFramebuffersEXT(1, &_FBO);

	_context = NULL;
	_validated = false;
      }

      /*! \brief Attaches this FBO as the current render target. */
      inline virtual void attach()
      {
	validate();
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _FBO);
	_context->setViewport(0, 0, getWidth(), getHeight());

	std::vector<GLenum> states(_colorTextures.size(), GL_NONE);
	for (size_t attachment(0); attachment < _colorTextures.size(); ++attachment)
	  if (_colorTextures[attachment])
	    states[attachment] = GL_COLOR_ATTACHMENT0_EXT + attachment;

	glDrawBuffers(states.size(), &states[0]);
      }

      /*! \brief Restores the screen FBO as the current render target. */
      inline 
      virtual void detach()
      {
	if (!_context)
	  M_throw() << "Cannot detach() an uninitialised FBO";

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);

	//Reset the draw states
	std::tr1::array<GLenum, 4> states 
	  = {{GL_FRONT_LEFT, GL_FRONT_RIGHT, GL_BACK_LEFT, GL_BACK_RIGHT}};
	glDrawBuffers(4, &states[0]);
      }

      inline void attachColorTexture(std::tr1::shared_ptr<Texture2D> coltex, size_t i)
      {
	if (!_context)
	  M_throw() << "Cannot attach textures to an uninitialised FBO";
	
	if (i >= _colorTextures.size())
	  M_throw() << "Out of range";

	_colorTextures[i] = coltex;
	_validated = false;
      }

      inline void attachDepthTexture(std::tr1::shared_ptr<Texture2D> depthtex)
      {
	if (!_context)
	  M_throw() << "Cannot attach textures to an uninitialised FBO";
	
	_depthTexture = depthtex;
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
	//First blit between the two FBO's
	glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, _FBO);
	glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER_EXT, other._FBO);
	glBlitFramebufferEXT(0, 0, getWidth(), getHeight(), 0, 0, 
			     other.getWidth(), other.getHeight(), opts, GL_NEAREST);
	
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
	  M_throw() << "Cannot fetch the color texture " << ID << " as the FBO has none bound";

	return *_colorTextures[ID];
      }

      /*! \brief Fetch the texture bound to the depth buffer. */
      inline Texture2D& getDepthTexture()
      { 
	if (!_depthTexture) 
	  M_throw() << "Cannot fetch the depth texture as the FBO has none bound";
	return *_depthTexture; 
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

	    //Bind the textures, or unbind the unbound textures ready for the completeness test
	    if (_depthTexture)
	      glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, 
					GL_TEXTURE_2D, _depthTexture->getGLHandle(), 0);
	    else
	      glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, 
					GL_TEXTURE_2D, 0, 0);

	    
	    std::vector<GLenum> states(_colorTextures.size(), GL_NONE);

	    for (size_t attachment(0); attachment < _colorTextures.size(); ++attachment)
	      if (_colorTextures[attachment])
		{
		  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT + attachment, GL_TEXTURE_2D, 
					    _colorTextures[attachment]->getGLHandle(), 0);
		  states[attachment] = GL_COLOR_ATTACHMENT0_EXT + attachment;
		}
	      else
		glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT + attachment, GL_TEXTURE_2D, 
					  0, 0);

	    glDrawBuffers(states.size(), &states[0]);

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
		M_throw() << "Failed to create FrameBufferObject: Unkown error code";
	      case GL_FRAMEBUFFER_COMPLETE_EXT: 
		break;
	      }
	    _validated = true;
	  }
      }

      Context* _context;
      std::vector<std::tr1::shared_ptr<Texture2D> > _colorTextures;
      std::tr1::shared_ptr<Texture2D> _depthTexture;
      GLuint _FBO;
      bool _validated;
    };
  }
}
