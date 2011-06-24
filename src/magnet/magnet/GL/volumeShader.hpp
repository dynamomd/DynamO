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
#include <magnet/math/vector.hpp>

namespace magnet {
  namespace GL {

    class volumeRenderer: public detail::Shader
    {
    public:      
      inline void build()
      {
	//First, call the build function in the shader
	Shader::build();

	_FocalLengthUniform = glGetUniformLocationARB(_shaderID,"FocalLength");
	_WindowSizeUniform = glGetUniformLocationARB(_shaderID,"WindowSize");
	_RayOriginUniform = glGetUniformLocationARB(_shaderID,"RayOrigin");
	_depthTexUniform = glGetUniformLocationARB(_shaderID,"DepthTexture");
	_nearUniform = glGetUniformLocationARB(_shaderID,"NearDist");
	_farUniform = glGetUniformLocationARB(_shaderID,"FarDist");
	_dataTexUniform = glGetUniformLocationARB(_shaderID,"DataTexture");
	_stepSizeUniform = glGetUniformLocationARB(_shaderID,"StepSize");
	_diffusiveLightingUniform = glGetUniformLocationARB(_shaderID,"DiffusiveLighting");
	_specularLightingUniform = glGetUniformLocationARB(_shaderID,"SpecularLighting");
	_ditherRayUniform = glGetUniformLocationARB(_shaderID,"DitherRay");
	_transferTexUniform = glGetUniformLocationARB(_shaderID,"TransferTexture");
      }

      inline void attach(GLfloat FocalLength, GLint width, GLint height, Vector Origin, 
			 GLint depthTex, GLint dataTex, GLint transferTex,
			 GLfloat NearDist, GLfloat FarDist,
			 GLfloat stepSize, GLfloat diff,
			 GLfloat spec, GLfloat dither)
      {
	glUseProgramObjectARB(_shaderID);
	glUniform1fARB(_FocalLengthUniform, FocalLength);
	glUniform2fARB(_WindowSizeUniform, width, height);
	glUniform3fARB(_RayOriginUniform, Origin[0], Origin[1], Origin[2]);
	glUniform1iARB(_depthTexUniform, depthTex);
	glUniform1fARB(_farUniform, FarDist);
	glUniform1fARB(_nearUniform, NearDist);
	glUniform1iARB(_dataTexUniform, dataTex);
	glUniform1fARB(_stepSizeUniform, stepSize);
	glUniform1fARB(_diffusiveLightingUniform, diff);
	glUniform1fARB(_specularLightingUniform, spec);
	glUniform1fARB(_ditherRayUniform, dither);
	glUniform1iARB(_transferTexUniform, transferTex);
      }

      virtual std::string vertexShaderSource();
      virtual std::string fragmentShaderSource();
      
    protected:
      GLuint _FocalLengthUniform;
      GLuint _WindowSizeUniform;
      GLuint _RayOriginUniform;
      GLuint _depthTexUniform;
      GLuint _nearUniform;
      GLuint _farUniform;
      GLuint _dataTexUniform;
      GLuint _stepSizeUniform;
      GLuint _diffusiveLightingUniform;
      GLuint _specularLightingUniform;
      GLuint _ditherRayUniform;
      GLuint _transferTexUniform;
    };
  }
}

#include <magnet/GL/detail/shaders/volumeShader.glh>
