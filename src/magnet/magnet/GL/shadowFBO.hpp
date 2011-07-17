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
    /*! \brief A Frame Buffer Object with appropriate depth texture
     * settings for a shadow mapping buffer.
     *
     */
    class shadowFBO : public FBO
    {
    public:
      virtual void init(GLsizei, GLsizei, GLint internalformat = GL_RGBA)
      { M_throw() << "Cannot use this initializer"; }

      /*! \brief Initializes the shadow FBO
       * 
       * \param length The side length of the FBO in pixels.
       */
      virtual void init(GLsizei length)
      {
	FBO::init(length, length);

	_depthTexture.bind(0);
	GLfloat l_ClampColor[] = {0.0, 0.0, 0.0, 0.0};
	glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, l_ClampColor);
	_depthTexture.parameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	_depthTexture.parameter(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	_depthTexture.parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
	_depthTexture.parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);

	//Enable shadow comparison
	_depthTexture.parameter(GL_TEXTURE_COMPARE_MODE, GL_COMPARE_R_TO_TEXTURE);
	//Shadow comparison should be true (ie not in shadow) if r<=texture
	_depthTexture.parameter(GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);
	//Shadow comparison should generate an INTENSITY result
	_depthTexture.parameter(GL_DEPTH_TEXTURE_MODE, GL_INTENSITY);
	
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _FBO);
	glDrawBuffer(GL_NONE);
	glReadBuffer(GL_NONE);
	// switch back to window-system-provided framebuffer
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
      }

      /*! \brief Sets up the FBO ready for the light perspective
       * render pass.
       */
      inline void setup()
      {
	//Use the fixed pipeline 
	glUseProgramObjectARB(0);
	//Render to the shadow maps FBO
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,_FBO);
	//Clear the depth buffer
	glClear(GL_DEPTH_BUFFER_BIT);
	//The viewport should change to the shadow maps size
	glViewport(0, 0, _width, _width);
      	//Draw back faces into the shadow map
	//glCullFace(GL_FRONT);	
	//Use flat shading for speed
	glShadeModel(GL_FLAT);
	//Mask color writes
	glColorMask(0, 0, 0, 0);
      }

      /*! \brief Restores the original screen FBO. */
      inline void restore()
      {
	//Restore the default FB
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);
	glShadeModel(GL_SMOOTH);
	glColorMask(1, 1, 1, 1);
      }
    };
  }
}
