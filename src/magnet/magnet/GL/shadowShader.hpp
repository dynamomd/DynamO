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

    class shadowShader: public detail::shader<shadowShader>
    {
    public:      
      inline void build()
      {
	//First, call the build function in the shader
	detail::shader<shadowShader>::build();
	
	//Now we fetch the uniforms out of the shader
	_shadowMapUniform 
	  = glGetUniformLocationARB(_shaderID,"ShadowMap");
	_shadowIntensityUniform 
	  = glGetUniformLocationARB(_shaderID,"shadowIntensity");
	_shadowMapStepXUniform 
	  = glGetUniformLocationARB(_shaderID,"xPixelOffset");
	_shadowMapStepYUniform 
	  = glGetUniformLocationARB(_shaderID,"yPixelOffset");
	_shadowMapEnable = glGetUniformLocationARB(_shaderID,"shadowMapping");
      }

      inline void attach(GLuint shadowTexture, size_t shadowSize, 
			 GLuint textureUnit, GLint shadowMapping, 
			 GLfloat intensity, size_t width, size_t height)
      {
	glUseProgramObjectARB(_shaderID);
	glUniform1iARB(_shadowMapUniform, textureUnit);
	glUniform1fARB(_shadowIntensityUniform, intensity);
	glUniform1fARB(_shadowMapStepXUniform, 1.0f / width);
	glUniform1fARB(_shadowMapStepYUniform, 1.0f / height);
	glUniform1iARB(_shadowMapEnable, shadowMapping);	
      }

      static inline std::string vertexShaderSource();
      static inline std::string fragmentShaderSource();
      
    protected:
      GLuint _shadowMapUniform;
      GLuint _shadowIntensityUniform;
      GLuint _shadowMapStepXUniform;
      GLuint _shadowMapStepYUniform;
      GLuint _shadowMapEnable;
    };
  }
}

#include <magnet/GL/detail/shaders/shadowShader.glh>
