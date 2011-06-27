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
#define STRINGIFY(A) #A

namespace magnet {
  namespace GL {
    
    class MultiplyTexture : public detail::Shader
    {
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

      virtual std::string vertexShaderSource()
      {
	return STRINGIFY( 
void main(void)
{
  gl_Position = ftransform();
  gl_TexCoord[0] = gl_MultiTexCoord0;
}
);
      }

      virtual std::string fragmentShaderSource()
      {
	return STRINGIFY(
uniform sampler2D u_Texture0; //input
uniform sampler2D u_Texture1; //Depth buffer

void main(void)
{
  gl_FragColor = texture2D(u_Texture0, gl_TexCoord[0].st) * texture2D(u_Texture1, gl_TexCoord[0].st);
}
);
      }
    };
  }
}

#undef STRINGIFY
