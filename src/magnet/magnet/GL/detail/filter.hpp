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

#include <string>
#include <magnet/GL/detail/shader.hpp>

namespace magnet {
  namespace GL {
    namespace detail {
      
      /* This is a CRTP base class that builds filters (convolution
       * kernels for textures).
       *
       * It requires that the type that inherits it, specifies its own
       * type in the template parameter (T) and defines static member
       * functions called  T::vertexShaderSource() and T::fragmentShaderSource().
       */
      template<class T>
      class filter<T> : public shader<filter<T> > 
      {
      public:
	//bind to an existing FBO
	filter(GLuint FBO, GLsizei width, GLsizei height):
	  _FBO(FBO),
	  _width(width),
	  _height(height)
	{}

	//Create our own FBO
	filter(GLsizei width, GLsizei height):
	  _width(width),
	  _height(height)
	{
	  glGenFramebuffersEXT(1, &_FBO);
	  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo->frame[i]);
	  	  
	  GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
	  if (status != GL_FRAMEBUFFER_COMPLETE_EXT)
	    M_throw() << "Error [" <<  status << "] while creating frame buffer";
	  
	  // switch back to screen framebuffer
	  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
	}
	
	
      protected:
	GLuint _FBO;
	
	GLsizei _width;
	GLsizei _height;
	
	//Helper functions to attach textures to the FBO
	void bindTexture(GLuint* texture)
	{
	  //Bind the FBO
	  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _FBO);

	  //Generate and bind the texture
	  glGenTextures(1, texture);
	  glBindTexture(GL_TEXTURE_2D, texture);
	  
	  //bind it to the FBO
	  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, fbo->texid[i], 0);
	  
	  //Set the textures properties
	  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F_ARB, _width, _height, 0, GL_RGBA, GL_FLOAT, NULL);
	  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	}

	void preInvoke()
	{
	  glUseProgram(_shaderID);
	  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _FBO);
	  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	  
	  glViewport(offsetX, offsetY, width, heigth);

	  //Save the matrix state
	  glMatrixMode(GL_PROJECTION);
	  glPushMatrix();
	  glLoadIdentity();

	  glMatrixMode(GL_MODELVIEW);
	  glPushMatrix();
	  glLoadIdentity();

	}

	//Before calling postInvoke, you should have called preInvoke,
	//bound the textures you wanted to, and set the shaders
	//uniforms.
	void postInvoke()
	{
	  //Draw a full screen quad to generate all the fragment shaders
	  glBegin(GL_QUADS);
	  glTexCoord2f(0.0f, 0.0f);
	  glVertex2d(-1, -1);
	  glTexCoord2f(1.0f, 0.0f);
	  glVertex2d(1, -1);
	  glTexCoord2f( 1.0f, 1.0f);
	  glVertex2d(1, 1);
	  glTexCoord2f(0.0f, 1.0f);
	  glVertex2d(-1, 1);
	  glEnd();

	  //Restore the matrix state
	  glMatrixMode(GL_PROJECTION);
	  glPopMatrix();

	  glMatrixMode(GL_MODELVIEW);
	  glPopMatrix();

	  //Restore the fixed pipeline
	  glUseProgramObjectARB(0);
	}

	void textureToScreen()
	{
	  //Save the matrix state and load identities
	  glMatrixMode(GL_PROJECTION);
	  glPushMatrix();
	  glLoadIdentity();

	  glMatrixMode(GL_MODELVIEW);
	  glPushMatrix();
	  glLoadIdentity();

	  
	  
	  //Restore the matrix state
	  glMatrixMode(GL_PROJECTION);
	  glPopMatrix();

	  glMatrixMode(GL_MODELVIEW);
	  glPopMatrix();
	}
      protected:
	
      };
    }
  }
}
