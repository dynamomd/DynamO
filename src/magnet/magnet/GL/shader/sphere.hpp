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
      /*! \brief A G-Buffer shader which billboards and raytraces spheres.
       */
      class SphereShader: public detail::Shader
      {
      public:
	virtual std::string initVertexShaderSource()
	{
	  return "#version 330\n" 
	    STRINGIFY(
uniform mat4 ViewMatrix;
uniform float global_scale;

layout (location = 0) in vec4 vPosition;
layout (location = 1) in vec4 vColor;
layout (location = 5) in vec4 iScale;

out vec4 color;
out float scale;

void main()
{
  color = vColor;
  scale = iScale.x * global_scale;
  gl_Position = ViewMatrix * vec4(vPosition.xyz, 1.0);
});
	}
	
	virtual std::string initGeometryShaderSource()
	{
	  return
	    "#version 330\n"
	    STRINGIFY(
uniform mat4 ProjectionMatrix;

layout(points) in;
layout(triangle_strip) out;
layout(max_vertices = 4) out;

in vec4 color[];
in float scale[];

flat out vec4 vert_color;
flat out float frag_scale;
flat out vec3 model_position_frag;
smooth out vec2 ordinate;

void VertexEmit(vec2 displacement, float amount)
{
  ordinate = displacement;
  vec4 shift = vec4(amount * displacement, 0.0, 0.0);
  vec4 eyespace_position = gl_in[0].gl_Position + shift;
  gl_Position = ProjectionMatrix * eyespace_position;
  EmitVertex();
}

void main()
{
  float halfsize = scale[0];

  //Standard data for each fragment
  vert_color = color[0];
  frag_scale = scale[0];
  model_position_frag = gl_in[0].gl_Position.xyz;
  VertexEmit(vec2(-1.0, -1.0), scale[0]);
  VertexEmit(vec2(-1.0, +1.0), scale[0]);
  VertexEmit(vec2(+1.0, -1.0), scale[0]);
  VertexEmit(vec2(+1.0, +1.0), scale[0]);
  EndPrimitive();
}
);
	}
	
	virtual std::string initFragmentShaderSource()
	{
	  return "#version 330\n"
	    STRINGIFY(
uniform mat4 ProjectionMatrix;

flat in vec4 vert_color;
flat in float frag_scale;
flat in vec3 model_position_frag;
smooth in vec2 ordinate;

layout (location = 0) out vec4 color_out;
layout (location = 1) out vec4 normal_out;
layout (location = 2) out vec4 position_out;

void main()
{
  //The ordinate variable contains the x and y position of the
  //sphere. Use the equation of a sphere to determine the z pos
  float z = 1.0 - ordinate.x * ordinate.x - ordinate.y * ordinate.y;

  //Discard the fragment if it lies outside the sphere
  if (z <= 0.0) discard;

  //Calculate the fragments real position on the sphere
  z = sqrt(z);
  vec4 frag_position_eye = vec4(model_position_frag + vec3(ordinate, z) * frag_scale, 1.0);
  //Calculate the fragments depth
  vec4 pos = ProjectionMatrix * frag_position_eye;
  gl_FragDepth = (pos.z / pos.w + 1.0) / 2.0; 

  //Write out the fragment's data
  position_out = frag_position_eye;
  color_out = vert_color;
  normal_out = vec4(ordinate.x, ordinate.y, z, 1.0);
});
	}
      };
    }
  }
}

#undef STRINGIFY
