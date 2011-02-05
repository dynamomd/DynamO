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
    
    class BilateralBlur : public detail::shader<BilateralBlur>
    {
    public:
      void build()
      {
	detail::shader<BilateralBlur>::build();

	glUseProgram(detail::shader<BilateralBlur>::_shaderID);

	_scaleUniform = glGetUniformLocationARB(detail::shader<BilateralBlur>::_shaderID,"scale");
	_totstrengthUniform = glGetUniformLocationARB(detail::shader<BilateralBlur>::_shaderID,"totStrength");

	_nearDistUniform = glGetUniformLocationARB(detail::shader<BilateralBlur>::_shaderID,"nearDist");
	_farDistUniform = glGetUniformLocationARB(detail::shader<BilateralBlur>::_shaderID,"farDist");

	_SSAOTextureUniform = glGetUniformLocationARB(detail::shader<BilateralBlur>::_shaderID,"u_Texture0");
	_depthTextureUniform = glGetUniformLocationARB(detail::shader<BilateralBlur>::_shaderID,"u_Texture2");

	glUseProgramObjectARB(0);
      }

      void invoke(GLint SSAOTextureID, GLint depthTextureID, GLuint _width, GLuint _height,
		  GLfloat pixelSkip, GLfloat totStrength, GLfloat neardist, GLfloat fardist)
      {
	//Setup the shader arguments
	glUseProgram(detail::shader<BilateralBlur>::_shaderID);
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

      static inline std::string vertexShaderSource();
      static inline std::string fragmentShaderSource();

    protected:
      GLint _SSAOTextureUniform, _depthTextureUniform;
      GLint _scaleUniform     ;
      GLint _totstrengthUniform;
      GLint _nearDistUniform, _farDistUniform;
    };
  }
}

#include <magnet/GL/detail/shaders/BilateralBlur.glh>
