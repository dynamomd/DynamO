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
      /*! \brief A deffered rendering (G-Buffer) shader which
	billboards/raytraces cylinders.

	This shader provides an extremely fast method to render
	perfect spheres in OpenGL. This method appears to outperform
	even the most poorly tesselated spheres. Only the position
	of the sphere (the input type is GL_POINTS) is needed as
	input (The radius of the sphere is passed in through the
	iScale vertex attribute). A geometry shader then converts
	each POINT into two triangles as a square bilboard. When the
	billboard is rasterized into fragments, in the fragment
	shader, each fragment is used to ray trace a sphere within
	the billboard. Thus, we only draw the front face of the
	sphere, using the absolute minimum input data, only two
	triangles at the cost of a slightly expensive fragment
	shader and an additional (trivial) geometry shader stage.

	Anti-aliasing can be achieved by forcing the GL state to
	evaluate all samples of the fragments using the
	GL_ARB_sample_shading extension when available. Something
	like \code glEnable(SAMPLE_SHADING_ARB);
	glMinSampleShadingARB(1.0); \endcode will enable
	multisampling on the spheres when possible.

	A discussion of this technique is given in the excellent
	online GL book by Jason L. McKesson at \url
	http://www.arcsynthesis.org/gltut/index.html in the chapter
	on lies and IMPOSTORS.
      */
      class CylinderShader: public detail::Shader
      {
      public:
	CylinderShader()
	{ 
	  defines("unshaded") = "false";
	}

	virtual std::string initVertexShaderSource()
	{
	  return STRINGIFY(
uniform mat4 ViewMatrix;
uniform float global_scale;

layout (location = 0) in vec4 vPosition;
layout (location = 1) in vec4 vColor;
layout (location = 4) in vec4 iOrientation;
layout (location = 5) in vec4 iScale;

out vec4 color;
out vec3 axis;
out float radius;
out float length;

vec3 qrot(vec4 q, vec3 v)
{ return v + 2.0 * cross(q.xyz, cross(q.xyz, v) + q.w * v); } 

void main()
{
  color = vColor;
  vec3 scale = iScale.xyz + vec3(equal(iScale.xyz, vec3(0.0,0.0,0.0)));
  radius = scale.x * global_scale * 0.5;
  length = scale.y * global_scale * 0.5;
  axis = normalize((ViewMatrix * vec4(qrot(iOrientation, vec3(0,0,1)), 0.0)).xyz);
  gl_Position = ViewMatrix * vec4(vPosition.xyz, 1.0);
});
	}
	
	virtual std::string initGeometryShaderSource()
	{
	  return STRINGIFY(
uniform mat4 ProjectionMatrix;

layout(points) in;
layout(triangle_strip) out;
layout(max_vertices = 4) out;

in vec4 color[];
in vec3 axis[];
in float radius[];
in float length[];

flat out vec4 vert_color;
flat out vec3 frag_axis;
flat out vec3 frag_center;
flat out vec2 frag_screen_perp;
flat out vec2 frag_screen_para;
smooth out vec2 frag_displacement;
flat out float frag_radius;
flat out float frag_length;

//Function to emit a bilboard vertex with all the correct output given
//the displacement
void VertexEmit(in vec2 displacement, in vec2 screen_perp, in vec2 screen_para)
{
  //The billboards need to be slightly larger to accommodate perspective changes.
  const float overdraw = 1.2;
  displacement *= overdraw;
  vec3 position = gl_in[0].gl_Position.xyz + length[0] * displacement.x * axis[0] + vec3(displacement.y * screen_perp + displacement.x * screen_para, 0.0);
  frag_axis = axis[0];
  frag_radius = radius[0];
  frag_length = length[0];
  vert_color = color[0];
  frag_center = gl_in[0].gl_Position.xyz;
  frag_screen_perp = screen_perp;
  frag_screen_para = screen_para;
  frag_displacement = displacement;
  gl_Position = ProjectionMatrix * vec4(position, gl_in[0].gl_Position.w);
  EmitVertex();
}

void main()
{
  //Standard data for each fragment
  float cosalpha = abs(dot(vec3(0.0,0.0,1.0), axis[0]));
  float da = radius[0] * cosalpha;
  float sinalpha = sqrt(1-cosalpha*cosalpha);
  vec2 screen_para = normalize(axis[0].xy);
  vec2 screen_perp = radius[0] * vec2(screen_para.y, -screen_para.x);
  screen_para *= da;
  VertexEmit(vec2(-1.0, -1.0), screen_perp, screen_para);
  VertexEmit(vec2(-1.0, +1.0), screen_perp, screen_para);
  VertexEmit(vec2(+1.0, -1.0), screen_perp, screen_para);
  VertexEmit(vec2(+1.0, +1.0), screen_perp, screen_para);
  EndPrimitive();
});
	}

	virtual std::string initFragmentShaderSource()
	{
	  return STRINGIFY(
uniform mat4 ProjectionMatrix;

flat in vec4 vert_color;
flat in vec3 frag_axis;
flat in vec3 frag_center;
flat in vec2 frag_screen_perp;
flat in vec2 frag_screen_para;
smooth in vec2 frag_displacement;
flat in float frag_radius;
flat in float frag_length;
layout (location = 0) out vec4 color_out;
layout (location = 1) out vec4 normal_out;
layout (location = 2) out vec4 position_out;

void main()
{
  vec3 fragment_position = frag_center + frag_length * frag_displacement.x * frag_axis + vec3(frag_displacement.y * frag_screen_perp + frag_displacement.x * frag_screen_para, 0.0);

)"\n#ifdef DRAWBILLBOARD\n"STRINGIFY(
  color_out = vert_color;
  normal_out = vec4(0.0);
  position_out = vec4(fragment_position, 1.0);
  vec4 pos = ProjectionMatrix * vec4(fragment_position, 1.0);
  gl_FragDepth = (pos.z / pos.w + 1.0) / 2.0;
)"\n#else\n"STRINGIFY(
   vec3 rij = frag_center;
   vec3 rij_planar = rij - dot(rij, frag_axis) * frag_axis;
   vec3 vij = normalize(-fragment_position);
   vec3 vij_planar = vij - dot(vij, frag_axis) * frag_axis;

   float A = dot(vij_planar, vij_planar);
   float B = dot(rij_planar, vij_planar);
   float C = dot(rij_planar, rij_planar) - frag_radius * frag_radius;
   float argument = B * B - A * C;
   if (argument < 0.0) discard;
   float t = - C / (B - sqrt(argument));
   vec3 hit = -t * vij;
   vec3 relative_hit = hit - frag_center;
   float axial_displacement = dot(relative_hit, frag_axis);
   if (abs(axial_displacement) > frag_length) discard;
   //normal_out = vec4(0.0); 
   normal_out = vec4(normalize(relative_hit - axial_displacement * frag_axis),1.0);
   color_out = vert_color;
   position_out = vec4(hit, 1.0);
   vec4 pos = ProjectionMatrix * vec4(hit, 1.0);
   gl_FragDepth = (pos.z / pos.w + 1.0) / 2.0;
)"\n#endif\n"STRINGIFY(
});
	}
      };

      /*! \brief A variant of the CylinderShader used for variance
          shadow mapping.
       */
      class CylinderVSMShader: public CylinderShader
      {
	virtual std::string initFragmentShaderSource()
	{
	  return STRINGIFY(
uniform mat4 ProjectionMatrix;

flat in float frag_radius;
flat in vec3 frag_center;
smooth in vec2 ordinate;

layout (location = 0) out vec4 color_out;

void main()
{
  vec3 billboard_frag_pos = frag_center + vec3(ordinate, 0.0) * frag_radius;
  vec3 ray_direction = normalize(billboard_frag_pos);

  float TD = dot(ray_direction, -frag_center);
  float c = dot(frag_center, frag_center) - frag_radius * frag_radius;
  float arg = TD * TD - c;
      
  if (arg < 0) discard;
  
  float t = - c / (TD - sqrt(arg));

  vec3 frag_position_eye = ray_direction * t;

  //Calculate the fragments depth
  vec4 pos = ProjectionMatrix * vec4(frag_position_eye, 1.0);
  gl_FragDepth = (pos.z / pos.w + 1.0) / 2.0;

  float depth = -frag_position_eye.z;
  float A = ProjectionMatrix[2].z;
  float B = ProjectionMatrix[3].z;
  float moment1 = 0.5 * (-A * depth + B) / depth + 0.5;
  float moment2 = moment1 * moment1;

  // Adjusting moments (this is sort of bias per pixel) using derivative
  float dx = dFdx(moment1);
  float dy = dFdy(moment1);
  moment2 += 0.25 * (dx * dx + dy * dy);
	
  color_out = vec4(moment1, moment2, 0.0, 1.0);
});
	}

      };
    }
  }
}

#undef STRINGIFY
