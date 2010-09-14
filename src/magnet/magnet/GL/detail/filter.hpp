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
      class filter : public shader<filter<T> > 
      {
      public:
	

	//bind to an existing FBO
	void build(GLuint FBO, GLsizei width, GLsizei height, GLint internalFormat, GLenum type)
	{
	  _FBO = FBO;
	  _width = width;
	  _height = height;
	  
	  bindTexture(internalFormat, type);
	}

	//Create our own FBO
	void build(GLsizei width, GLsizei height, GLint internalFormat, GLenum type)
	{
	  _width = width;
	  _height = height;

	  glGenFramebuffersEXT(1, &_FBO);
	  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _FBO);

	  bindTexture(internalFormat, type);
	  	  
	  // switch back to screen framebuffer
	  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
	}
	
	void renderOutput()
	{
	  //Save the matrix state
	  glMatrixMode(GL_PROJECTION);
	  glPushMatrix();
	  glLoadIdentity();

	  glMatrixMode(GL_MODELVIEW);
	  glPushMatrix();
	  glLoadIdentity();

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
	}
	
      protected:
	GLuint _FBO, _outputTexture;
	GLsizei _width;
	GLsizei _height;
	
	void bindTexture(GLint internalFormat, GLenum type)
	{
	  glGenTextures(1, &_outputTexture);	
	  glBindTexture(GL_TEXTURE_2D, _outputTexture);
	  
	  glTexImage2D(GL_TEXTURE_2D, 0, internalFormat, _width, _height, 0, GL_RGBA, type, NULL);
	  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, _outputTexture, 0);
	  
	  GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
	  if (status != GL_FRAMEBUFFER_COMPLETE_EXT)
	    M_throw() << "Error [" <<  status << "] while creating frame buffer";
	}

	void preInvoke()
	{
	  glUseProgram(shader<filter<T> >::_shaderID);
	  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _FBO);
	  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	  
	  glViewport(0, 0, _width, _height);

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

      protected:
	
      };
    }
  }
}
