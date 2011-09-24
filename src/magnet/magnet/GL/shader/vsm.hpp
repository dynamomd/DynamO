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
      /*! \brief An instancing depth only shader for generating
	variance shadow map textures.
	
	For more information on variance shadow mapping, consult the
	original paper http://www.punkuser.net/vsm/
       */
      class VSMShader: public detail::Shader
      {
      public:
	virtual std::string initVertexShaderSource()
	{
	  return "#version 330\n"
	    STRINGIFY(
uniform mat4 ProjectionMatrix;
uniform mat4 ViewMatrix;

layout (location = 0) in vec4 vPosition;
layout (location = 3) in vec4 iOrigin;
layout (location = 4) in vec4 iOrientation;
layout (location = 5) in vec4 iScale;

smooth out float depth;

vec3 qrot(vec4 q, vec3 v)
{ return v + 2.0 * cross(cross(v,q.xyz) + q.w * v, q.xyz); }

void main()
{
  //Rotate the vertex according to the instance transformation, and
  //then move it to the instance origin.
  vec3 scale = iScale.xyz + vec3(equal(iScale.xyz, vec3(0.0))) * iScale.x;
  vec4 vVertex = ViewMatrix * vec4(qrot(iOrientation, vPosition.xyz * scale)
				   + iOrigin.xyz, 1.0);
  depth = -vVertex.z;
  vec4 pos = ProjectionMatrix * vVertex;
  gl_Position = pos;
});
	}
	
	virtual std::string initFragmentShaderSource()
	{ return "#version 330\n"
	    STRINGIFY(
uniform mat4 ProjectionMatrix;
smooth in float depth;

layout (location = 0) out vec4 color_out;
void main()
{
  float A = ProjectionMatrix[2].z;
  float B = ProjectionMatrix[3].z;
  float moment1 = 0.5 * (-A * depth + B) / depth + 0.5;
  float moment2 = moment1 * moment1;

  // Adjusting moments (this is sort of bias per pixel) using derivative
  float dx = dFdx(moment1);
  float dy = dFdy(moment1);
  moment2 += 0.25 * (dx * dx + dy * dy);
	
  color_out = vec4(moment1, moment2, 0, 1.0);
});
	}
      };
    }
  }
}

#undef STRINGIFY
