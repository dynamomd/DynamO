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
    class shadowFBO : public FBO
    {
    public:
      virtual void init(GLsizei, GLsizei, GLint, GLenum, GLenum)
      { M_throw() << "Cannot use this initializer"; }

      virtual void init(GLsizei length)
      {
	FBO::init(length, length);
	glBindTexture(GL_TEXTURE_2D, _depthTexture);

	GLfloat l_ClampColor[] = {0.0, 0.0, 0.0, 0.0};
	glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, l_ClampColor);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);

	//Enable shadow comparison
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_COMPARE_R_TO_TEXTURE);
	//Shadow comparison should be true (ie not in shadow) if r<=texture
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);
	//Shadow comparison should generate an INTENSITY result
	glTexParameteri(GL_TEXTURE_2D, GL_DEPTH_TEXTURE_MODE, GL_INTENSITY);
	
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _FBO);
	glDrawBuffer(GL_NONE);
	glReadBuffer(GL_NONE);
	// switch back to window-system-provided framebuffer
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
      }
      
      inline 
      virtual void resize(GLsizei width, GLsizei height)
      {
	if (width  != height) 
	  M_throw() << "Shadow maps should be square!";

	//Skip identity operations
	if (width  == _width) return;

	FBO::resize(width, height);
      }

      inline void setup()
      {
	//Use the fixed pipeline 
	glUseProgramObjectARB(0);

	//Render to the shadow maps FBO
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,_FBO);
	//Clear the depth buffer
	glClear(GL_DEPTH_BUFFER_BIT);

	//Save state
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	
	//The viewport should change to the shadow maps size
	glViewport(0, 0, _width, _width);
      	//Draw back faces into the shadow map
	//glCullFace(GL_FRONT);	
	//Use flat shading for speed
	glShadeModel(GL_FLAT);
	//Mask color writes
	glColorMask(0, 0, 0, 0);
      }

      inline void restore()
      {
	//Restore the draw mode
	glPopAttrib();

	//Restore the default FB
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);
      }
    };
  }
}
