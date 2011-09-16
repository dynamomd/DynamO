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

attribute vec4 vPosition;
attribute vec4 vColor;
attribute vec4 vNormal;
attribute vec4 iOrigin;
attribute vec4 iOrientation;
attribute vec4 iScale;

varying vec4 color;
varying vec3 normal;
varying vec3 position;

vec3 qrot(vec4 q, vec3 v)
{ return v + 2.0 * cross(cross(v,q.xyz) + q.w * v, q.xyz); } 

void main()
{
  color = vColor;
  //We store the normals in object space
  normal = normalize(qrot(iOrientation, vNormal.xyz));
  
  vec3 scale = iScale.xyz + vec3(equal(iScale.xyz, vec3(0.0))) * iScale.x;
  vec3 vVertex = qrot(iOrientation, vPosition.xyz * scale) + iOrigin.xyz;

  //We store the model space position of the vertex
  position = vVertex;

  gl_Position = ProjectionMatrix * (ViewMatrix * vec4(vVertex, 1.0));
});
	}
	
	virtual std::string initFragmentShaderSource()
	{
	  return STRINGIFY(
varying vec4 color;
varying vec3 normal;
varying vec3 position;

void main()
{
  gl_FragData[0] = vec4(color.rgb, 1.0);

  vec3 outnormal = normal;
  if (!gl_FrontFacing) outnormal = -outnormal;
  gl_FragData[1] = vec4(outnormal, 1.0);

  gl_FragData[2] = vec4(position.xyz, 1.0);
});
	}
      };
    }
  }
}

#undef STRINGIFY
