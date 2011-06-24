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
#include <algorithm>

namespace magnet {
  namespace GL {
    
    class SSAO : public detail::Shader
    {
    public:
      void build()
      {
	detail::Shader::build();

	glUseProgram(_shaderID);

	_radiusUniform      = glGetUniformLocationARB(_shaderID,"radius");
	_totstrengthUniform = glGetUniformLocationARB(_shaderID,"totStrength");
	_depthDropoffUniform    = glGetUniformLocationARB(_shaderID,"depthDropoff");
	_offsetUniform      = glGetUniformLocationARB(_shaderID,"offset");

	_nearDistUniform = glGetUniformLocationARB(_shaderID,"nearDist");
	_farDistUniform = glGetUniformLocationARB(_shaderID,"farDist");

	_colorTextureUniform = glGetUniformLocationARB(_shaderID,"u_Texture0");
	_normalTextureUniform = glGetUniformLocationARB(_shaderID,"u_Texture1");
	_depthTextureUniform = glGetUniformLocationARB(_shaderID,"u_Texture2");
	_rnmTextureUniform = glGetUniformLocationARB(_shaderID, "rnm");

	glUseProgramObjectARB(0);
      }

      void invoke(GLint colorTextureID, GLint normalTextureID, 
		  GLint depthTextureID, GLint rnmTextureID, 
		  GLuint _width, GLuint _height, 
		  GLfloat radius, GLfloat totStrength, GLfloat strength, 
		  size_t randomTexSize, GLfloat neardist, GLfloat fardist)
      {
	//Setup the shader arguments
	glUseProgram(_shaderID);
	//Horizontal application
	glUniform1iARB(_colorTextureUniform, colorTextureID);
	glUniform1iARB(_normalTextureUniform, normalTextureID);
	glUniform1iARB(_depthTextureUniform, depthTextureID);
	glUniform1iARB(_rnmTextureUniform, rnmTextureID);

	glUniform1fARB(_radiusUniform, radius);
	glUniform1fARB(_totstrengthUniform, totStrength);
	glUniform1fARB(_depthDropoffUniform, strength);

	glUniform1fARB(_nearDistUniform, neardist);
	glUniform1fARB(_farDistUniform, fardist);

	//Just tile the texture over the screen by multiplying the texture reads by this factor
	GLfloat offset = GLfloat(std::max(_width, _height)) / randomTexSize; 
	glUniform1fARB(_offsetUniform, offset);
	
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
      GLint _colorTextureUniform, _normalTextureUniform, _depthTextureUniform, _rnmTextureUniform;
      GLint _nearDistUniform, _farDistUniform;
      GLint _radiusUniform;
      GLint _totstrengthUniform;
      GLint _depthDropoffUniform;
      GLint _offsetUniform;
    };
  }
}

#include <magnet/GL/detail/shaders/SSAO.glh>
