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
    class shadowFBO
    {
    public:
      inline shadowFBO():
	_length(0)
      {}

      inline void init(GLsizei length)
      {
	if (_length) 
	  M_throw() << "shadowFBO has already been initialised!";

	if (!length)
	  M_throw() << "Cannot initialise a shadowFBO with a side length == 0!";

	_length = length;

	glGenFramebuffersEXT(1, &_FBO);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _FBO);
		
	//bind the depth texture
	glGenTextures(1, &_depthTexture);
	glBindTexture(GL_TEXTURE_2D, _depthTexture);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT24, _length, _length, 0,
		     GL_DEPTH_COMPONENT, GL_UNSIGNED_INT, NULL);
	
	GLfloat l_ClampColor[] = {1.0, 1.0,1.0,1.0};
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
	
	// attach the texture to FBO depth attachment point
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, 
				  GL_TEXTURE_2D, _depthTexture, 0);
	glDrawBuffer(GL_NONE);
	glReadBuffer(GL_NONE);

	// check FBO status
	GLenum FBOstatus = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
	if(FBOstatus != GL_FRAMEBUFFER_COMPLETE_EXT)
	  M_throw() << "GL_FRAMEBUFFER_COMPLETE_EXT failed";
	
	// switch back to window-system-provided framebuffer
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
      }
      
      inline ~shadowFBO()
      {
	if (_length)
	  {
	    glDeleteFramebuffersEXT(1, &_FBO);
	    glDeleteTextures(1, &_depthTexture);
	  }
      }

      inline void setup(const lightInfo& light)
      {
	//Now render the scene from the lights perspective

	//Use the fixed pipeline 
	glUseProgramObjectARB(0);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadMatrixf(light._projectionMatrix);
	
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadMatrixf(light._viewMatrix);
	

	//Render to the shadow maps FBO
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,_FBO);
	//Clear the depth buffer
	glClear(GL_DEPTH_BUFFER_BIT);

	//Save state
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	
	//The viewport should change to the shadow maps size
	glViewport(0, 0, _length, _length);
      	//Draw back faces into the shadow map
	glCullFace(GL_FRONT);	
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
	
	////Restore the Camera matricies
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();	
      }

      inline GLuint getFBO() { return _FBO; }
      inline GLuint getShadowTexture() { return _depthTexture; }
      inline GLsizei getLength() { return _length; }
    private:
      GLuint _FBO, _depthTexture;
      GLsizei _length;

    };
  }
}
