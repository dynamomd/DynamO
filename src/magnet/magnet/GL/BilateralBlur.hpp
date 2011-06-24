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

#include <magnet/GL/detail/shader.hpp>

namespace magnet {
  namespace GL {
    
    class BilateralBlur : public detail::Shader
    {
    public:
      void build()
      {
	detail::Shader::build();

	glUseProgram(_shaderID);
	_scaleUniform = glGetUniformLocationARB(_shaderID,"scale");
	_totstrengthUniform = glGetUniformLocationARB(_shaderID,"totStrength");

	_nearDistUniform = glGetUniformLocationARB(_shaderID,"nearDist");
	_farDistUniform = glGetUniformLocationARB(_shaderID,"farDist");

	_SSAOTextureUniform = glGetUniformLocationARB(_shaderID,"u_Texture0");
	_depthTextureUniform = glGetUniformLocationARB(_shaderID,"u_Texture2");

	glUseProgramObjectARB(0);
      }

      void invoke(GLint SSAOTextureID, GLint depthTextureID, GLuint _width, GLuint _height,
		  GLfloat pixelSkip, GLfloat totStrength, GLfloat neardist, GLfloat fardist)
      {
	//Setup the shader arguments
	glUseProgram(_shaderID);
	//Horizontal application
	glUniform1iARB(_SSAOTextureUniform, SSAOTextureID);
	glUniform1iARB(_depthTextureUniform, depthTextureID);

	glUniform2fARB(_scaleUniform, pixelSkip / _width, pixelSkip / _height);
	glUniform1fARB(_totstrengthUniform, totStrength);

	glUniform1fARB(_nearDistUniform, neardist);
	glUniform1fARB(_farDistUniform, fardist);
	  
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

      virtual std::string vertexShaderSource();
      virtual std::string fragmentShaderSource();

    protected:
      GLint _SSAOTextureUniform, _depthTextureUniform;
      GLint _scaleUniform     ;
      GLint _totstrengthUniform;
      GLint _nearDistUniform, _farDistUniform;
    };
  }
}

#include <magnet/GL/detail/shaders/BilateralBlur.glh>
