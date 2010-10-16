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

#include <magnet/GL/detail/shader.hpp>

namespace magnet {
  namespace GL {
    
    class SSAO : public detail::shader<SSAO>
    {
    public:
      void build()
      {
	detail::shader<SSAO>::build();

	//Get the shader args
	glUseProgram(detail::shader<SSAO>::_shaderID);
	_radiusUniform = glGetUniformLocationARB(detail::shader<SSAO>::_shaderID,"radius");
	_colorTextureUniform = glGetUniformLocationARB(detail::shader<SSAO>::_shaderID,"u_Texture0");
	_depthTextureUniform = glGetUniformLocationARB(detail::shader<SSAO>::_shaderID,"u_Texture1");
	glUseProgramObjectARB(0);
      }

      void invoke(GLint colorTextureID, GLint depthTextureID, GLuint _width, GLuint _height, GLfloat radius)
      {
	//Setup the shader arguments
	glUseProgram(detail::shader<SSAO>::_shaderID);
	//Horizontal application
	glUniform1iARB(_colorTextureUniform, colorTextureID);
	glUniform1iARB(_depthTextureUniform, depthTextureID);

	glUniform1fARB(_radiusUniform, radius);
	  
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
      GLint _colorTextureUniform, _depthTextureUniform;
      GLint _radiusUniform;
    };
  }
}

#include <magnet/GL/detail/shaders/SSAO.glh>
