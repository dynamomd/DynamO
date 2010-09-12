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

    class shadowShader: public detail::shader<shadowShader>
    {
    public:      
      
      inline void build()
      {
	//First, call the build function in the shader
	detail::shader<shadowShader>::build();
	
	//Now we fetch the uniforms out of the shader
	_shadowMapUniform = glGetUniformLocationARB(_shaderID,"ShadowMap");
	_shadowMapStepXUniform = glGetUniformLocationARB(_shaderID,"xPixelOffset");
	_shadowMapStepYUniform = glGetUniformLocationARB(_shaderID,"yPixelOffset");
      }

      inline void attach(GLuint shadowTexture, size_t shadowSize)
      {
	//Apply a matrix to the texture 7 unit
	double modelView[16];
	double projection[16];
	
	// This is matrix transform every coordinate x,y,z
	// x = x* 0.5 + 0.5 
	// y = y* 0.5 + 0.5 
	// z = z* 0.5 + 0.5 
	// Moving from unit cube [-1,1] to [0,1]  
	const GLdouble bias[16] = {	
	  0.5, 0.0, 0.0, 0.0, 
	  0.0, 0.5, 0.0, 0.0,
	  0.0, 0.0, 0.5, 0.0,
	  0.5, 0.5, 0.5, 1.0};
	
	// Grab modelview and transformation matrices
	glGetDoublev(GL_MODELVIEW_MATRIX, modelView);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);
		
	glMatrixMode(GL_TEXTURE);
	glActiveTextureARB(GL_TEXTURE7);
	
	glLoadIdentity();	
	glLoadMatrixd(bias);
	
	// concatating all matrice into one.
	glMultMatrixd (projection);
	glMultMatrixd (modelView);
	
	// Go back to normal matrix mode
	glMatrixMode(GL_MODELVIEW);

	glUseProgramObjectARB(_shaderID);
	glUniform1iARB(_shadowMapUniform,7);
	glUniform1fARB(_shadowMapStepXUniform, 1.0 / (shadowSize * 2));
	glUniform1fARB(_shadowMapStepYUniform, 1.0 / (shadowSize * 2));
	glActiveTextureARB(GL_TEXTURE7);
	glBindTexture(GL_TEXTURE_2D, shadowTexture);
      }

      GLuint _shadowMapUniform;
      GLuint _shadowMapStepXUniform;
      GLuint _shadowMapStepYUniform;

      static inline std::string vertexShaderSource();
      static inline std::string fragmentShaderSource();
    };
  }
}

#include <magnet/GL/detail/shaders/shadowShader.glh>
