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
      class DepthRenderShader: public detail::Shader
      {
      public:
	DepthRenderShader(): Shader(true, true) {}

	virtual std::string initVertexShaderSource()
	{
	  return STRINGIFY(
uniform mat4 ProjectionMatrix;
uniform mat4 ViewMatrix;

attribute vec4 vPosition;
attribute vec4 iOrigin;
attribute vec4 iOrientation;
attribute vec4 iScale;

////Quaternion mathmatics
//https://mollyrocket.com/forums/viewtopic.php?p=6154
vec3 qrot(vec4 q, vec3 v)
{
  return v + 2.0 * cross(cross(v,q.xyz) + q.w * v, q.xyz);
}

void main()
{
  //Rotate the vertex according to the instance transformation, and
  //then move it to the instance origin.
  vec4 vVertex = ViewMatrix * vec4(qrot(iOrientation, vPosition.xyz * iScale.xyz)
				   + iOrigin.xyz, 1.0);
  gl_Position = ProjectionMatrix * vVertex;
});
	}
	
	virtual std::string initFragmentShaderSource()
	{ return STRINGIFY(void main() {}); }
      };
    }
  }
}

#undef STRINGIFY
