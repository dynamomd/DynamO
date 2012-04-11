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
#include <iostream>
#define STRINGIFY(A) #A

namespace magnet {
  namespace GL {
    namespace shader {
      /*! \brief An instancing depth only shader for generating shadow maps.
       */
      class SimpleRenderShader: public detail::Shader
      {
      public:
	virtual std::string initVertexShaderSource()
	{
	  return STRINGIFY( 
uniform mat4 ProjectionMatrix;
uniform mat4 ViewMatrix;

layout (location = 0) in vec4 vPosition;
layout (location = 1) in vec4 vColor;
layout (location = 3) in vec4 iOrigin;
layout (location = 4) in vec4 iOrientation;
layout (location = 5) in vec4 iScale;

flat out vec4 color;

vec3 qrot(vec4 q, vec3 v)
{ return v + 2.0 * cross(cross(v,q.xyz) + q.w * v, q.xyz); } 

void main()
{
  color = vColor;
  
  vec3 scale = iScale.xyz + vec3(equal(iScale.xyz, vec3(0.0))) * iScale.x;
  vec4 vVertex
    = ViewMatrix
    * vec4(qrot(iOrientation, vPosition.xyz * scale) + iOrigin.xyz, 1.0);

  gl_Position = ProjectionMatrix * vVertex;
});
	}
	
	virtual std::string initFragmentShaderSource()
	{ return STRINGIFY(
flat in vec4 color;
layout (location = 0) out vec4 color_out;
void main()
{ color_out = color; }
); 
	}
      };
    }
  }
}

#undef STRINGIFY
