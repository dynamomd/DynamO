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
#include <magnet/GL/shader/detail/shader.hpp>
#define STRINGIFY(A) #A

namespace magnet {
  namespace GL {
    namespace shader {
      class NormalShader : public detail::Shader
      {
      public:
	inline void attach()
	{
	  //Setup the shader arguments
	  glUseProgram(_shaderID);
	  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	}

	virtual std::string vertexShaderSource()
	{
	  return STRINGIFY(
varying vec3 Normal;
void main()
{
  vec4 viewPos = gl_ModelViewMatrix * gl_Vertex;
  gl_Position = ftransform();
  Normal = normalize((gl_ModelViewMatrix * vec4(gl_Normal.xyz,0.0)).xyz);
});
	}
	
	virtual std::string fragmentShaderSource()
	{
	  return STRINGIFY(
varying vec3 Normal;
void main( void )
{
  gl_FragColor = vec4(0.5 * normalize(Normal) + 0.5, 1.0);
});
	}
      };
    }
  }
}

#undef STRINGIFY
