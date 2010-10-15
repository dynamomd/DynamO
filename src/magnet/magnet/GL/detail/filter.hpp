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
       * It requires that the type that inherits it, specifies the
       * width of the seperable filter (stencilwidth) as a template
       * and a suitable row (2*stencilwidth) of floats of the filter
       * block in a static const GLfloat ::weights member variable.
       */
      template<class T, int stencilwidth>
      class filter : public shader<filter<T,stencilwidth> >
      {
      public:

	void build()
	{
	  shader<filter<T, stencilwidth> >::build();

	  //Get the shader args
	  glUseProgram(shader<filter<T, stencilwidth> >::_shaderID);
	  _scaleUniform = glGetUniformLocationARB(shader<filter<T, stencilwidth> >::_shaderID,"u_Scale");	  
	  _textureUniform = glGetUniformLocationARB(shader<filter<T, stencilwidth> >::_shaderID,"u_Texture0");	  

	  //Set the weights now
	  GLint weightsUniform = glGetUniformLocationARB(shader<filter<T, stencilwidth> >::_shaderID, "weights");
	  glUniform1fvARB(weightsUniform, stencilwidth * stencilwidth, T::weights());

	  //Restore the fixed pipeline
	  glUseProgramObjectARB(0);
	}

	void invoke(GLint TextureID, GLuint _width, GLuint _height)
	{
	  
	  //Setup the shader arguments
	  glUseProgram(shader<filter<T, stencilwidth> >::_shaderID);
	  //Horizontal application
	  glUniform2fARB(_scaleUniform, 1.0 / _width, 1.0 / _height);
	  glUniform1iARB(_textureUniform, TextureID);
	  
	  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	  //Set the viewport
	  glPushAttrib(GL_VIEWPORT_BIT);
	  glViewport(0, 0, _width, _height);

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

	  //Restore the viewport
	  glPopAttrib();

	  //Restore the fixed pipeline
	  glUseProgramObjectARB(0);
	}
	
	static inline std::string vertexShaderSource();
	static inline std::string fragmentShaderSource();

      protected:
	GLint _scaleUniform, _textureUniform;
      };
    }
  }
}

#include <magnet/GL/detail/shaders/filter.glh>
