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

#include <string>
#include <magnet/GL/detail/shader.hpp>

namespace magnet {
  namespace GL {
    namespace detail {
      
      class SSFilter : public Shader
      {
      protected:
	void build(int stencilwidth)
	{
	  _stencilwidth = stencilwidth;

	  Shader::build();
	  //Get the shader args
	  glUseProgram(_shaderID);
	  _scaleUniform = glGetUniformLocationARB(_shaderID,"u_Scale");	  
	  _textureUniform = glGetUniformLocationARB(_shaderID,"u_Texture0");	  

	  //Set the weights now
	  GLint weightsUniform = glGetUniformLocationARB(_shaderID, "weights");
	  glUniform1fvARB(weightsUniform, stencilwidth * stencilwidth, weights());

	  //Restore the fixed pipeline
	  glUseProgramObjectARB(0);
	}

      public:
	void invoke()
	{
	  //Setup the shader arguments
	  glUseProgram(_shaderID);	  
	  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	  drawScreenQuad();
	  //Restore the fixed pipeline
	  glUseProgramObjectARB(0);
	}
	
	virtual std::string vertexShaderSource();
	virtual std::string fragmentShaderSource();

      protected:
	virtual const GLfloat* weights() = 0;

	GLint _scaleUniform, _textureUniform;
	int _stencilwidth;
      };
    }
  }
}

#include <magnet/GL/detail/shaders/filter.glh>
