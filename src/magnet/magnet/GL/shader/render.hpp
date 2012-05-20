/*    dynamo:- Event driven molecular dynamics simulator 
 *    http://www.dynamomd.org
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
      /*! \brief A G-Buffer render shader.

	This shader outputs all of the information needed for deffered
	shading calculations later on.
       */
      class RenderShader: public detail::Shader
      {
      public:
	virtual std::string initVertexShaderSource()
	{
	  return STRINGIFY(
uniform mat4 ProjectionMatrix;
uniform mat4 ViewMatrix;

layout (location = 0) in vec4 vPosition;
layout (location = 1) in vec4 vColor;
layout (location = 2) in vec4 vNormal;
layout (location = 3) in vec4 iOrigin;
layout (location = 4) in vec4 iOrientation;
layout (location = 5) in vec4 iScale;

flat out vec4 color;
smooth out vec3 normal;
smooth out vec3 position;

vec3 qrot(vec4 q, vec3 v)
{ return v + 2.0 * cross(cross(v,q.xyz) + q.w * v, q.xyz); } 

void main()
{
  color = vColor;
  //We store the normals in eye-space. The w coordinate is 0 to
  //prevent translations having any effect. The ViewMatrix must have
  //no scaling, only translations and rotations.
  normal = (ViewMatrix * vec4(qrot(iOrientation, vNormal.xyz), 0.0)).xyz;
  
  vec3 scale = iScale.xyz + vec3(equal(iScale.xyz, vec3(0.0))) * iScale.x;
  vec4 vVertex = ViewMatrix
    * vec4(qrot(iOrientation, vPosition.xyz * scale) + iOrigin.xyz, 1.0);

  //We store the eye-space position of the vertex
  position = vVertex.xyz;
  gl_Position = ProjectionMatrix * vVertex;
});
	}
	
	virtual std::string initFragmentShaderSource()
	{
	  return STRINGIFY(
flat in vec4 color;
smooth in vec3 normal;
smooth in vec3 position;

layout (location = 0) out vec4 color_out;
layout (location = 1) out vec4 normal_out;
layout (location = 2) out vec4 position_out;

void main()
{
  color_out = color;

  vec3 outnormal = (!gl_FrontFacing) ? -normal : normal;

  float nrm_length = length(outnormal);
  nrm_length += float(nrm_length == 0);
  outnormal /= nrm_length;

  normal_out = vec4(outnormal, 1.0);
  position_out = vec4(position, 1.0);
});
	}
      };
    }
  }
}

#undef STRINGIFY
