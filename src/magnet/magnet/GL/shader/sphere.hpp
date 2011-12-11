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
      /*! \brief A deffered rendering (G-Buffer) shader which
          billboards/raytraces spheres.

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

	  A final note is that we put the billboards exactly one
	  radius towards the viewer. This means the fragment shader
	  only increases the fragments depth after raytracing. This
	  coupled with the use of the GL_ARB_conservative_depth
	  extension allows early-Z culling to still be performed for
	  the standard GL_LESS/GL_LEQUAL depth testing!

	  There's tricky trick here too. The billboard is sized and
	  passed through the Projection matrix at the center depth of
	  the ball. This is all well and good, we will get a billboard
	  at exactly the correct xy size on the screen. But then we
	  wish to set the depth to exactly +r in projected space so
	  that the billboard depth is at the minimum the sphere might
	  yield. We can't just project the billboard at (position.z+r)
	  in the first place as with a perspecitve, the billboard will
	  be larger than requried to render the sphere! We have to
	  calculate the screen space z we want, project the billboard
	  in the middle and then replace its z component with the
	  screen space z we want multiplied by the w coordinate! See
	  the code for more details.

	  The only problem with this whole approach is that the
	  spheres aren't "fisheyed" correctly when they are very close
	  to the screen, when the triangular tesselated spheres are
	  fisheyed. I can't fix this without doing a better raytrace
	  in the fragment shader, which I haven't bothered to figure
	  out.
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

//Function to emit a bilboard vertex with all the correct output given
//the displacement
void VertexEmit(in vec2 displacement, in float finalz)
{
  ordinate = displacement;

  vec4 shift = vec4(scale[0] * displacement, 0.0, 0.0);
  vec4 eyespace_position = gl_in[0].gl_Position + shift;
  vec4 proj_position =  ProjectionMatrix * eyespace_position;

  //Here we do the cryptic thing of inputting the finalz times by the
  //w coordinate, so that when the w divide is performed by OpenGL, we
  //will recover the w divide.
  gl_Position = vec4(proj_position.xy, finalz * proj_position.w, proj_position.w);
  EmitVertex();
}

void main()
{
  //All of this paragraph concerns the manipulations required to get
  //the GL_ARB_conservative_depth extension to work without increasing
  //the render load.

  //We want to place the billboard at the front of the sphere, but we
  //want to make its screen dimensions those of a billboard halfway
  //through the sphere. Unfortunately, if we have perspective,
  //altering the depth will alter the screen dimensions too.

  //To get around this, we must calculate the
  //post-perspective/w-divided z we want (called finalz). We then
  //place our billboard vertices at the sphere middle but replace its
  //z component with finalz multiplied by the sphere middle w value,
  //so when OpenGL does the w-divide, we recover the correct depth!
  //Horray!

  //The projected position of the center of the front of the sphere
  vec4 centerpos = ProjectionMatrix * (gl_in[0].gl_Position + vec4(0.0, 0.0, scale[0], 0.0));
  
  //We now perform the w divide so that we have the minimum
  //screen-space depth of the sphere.
  float finalz = centerpos.z / centerpos.w;

  //Standard data for each fragment
  vert_color = color[0];
  frag_scale = scale[0];
  model_position_frag = gl_in[0].gl_Position.xyz;
  VertexEmit(vec2(-1.0, -1.0), finalz);
  VertexEmit(vec2(-1.0, +1.0), finalz);
  VertexEmit(vec2(+1.0, -1.0), finalz);
  VertexEmit(vec2(+1.0, +1.0), finalz);
  EndPrimitive();
}
);
	}

	virtual std::string initFragmentShaderSource()
	{
	  return "#version 330\n"
	    "#ifdef GL_ARB_conservative_depth\n"
	    "#extension GL_ARB_conservative_depth : enable\n"
	    "layout (depth_greater) out float gl_FragDepth;"
	    "#endif\n"
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
  float z = 1.0 - dot(ordinate,ordinate);

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

      /*! \brief A variant of the SphereShader used for variance
          shadow mapping.
       */
      class SphereVSMShader: public detail::Shader
      {
	virtual std::string initVertexShaderSource()
	{
	  return "#version 330\n" 
	    STRINGIFY(
uniform mat4 ViewMatrix;
uniform float global_scale;

layout (location = 0) in vec4 vPosition;
layout (location = 5) in vec4 iScale;
out float scale;
out float depth;

void main()
{
  scale = iScale.x * global_scale;
  vec4 vVertex = ViewMatrix * vec4(vPosition.xyz, 1.0); 
  gl_Position = vVertex;
  depth = -vVertex.z;
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

in float scale[];
in float depth[];

flat out float frag_scale;
flat out float frag_depth;
flat out vec3 model_position_frag;
smooth out vec2 ordinate;

//Function to emit a bilboard vertex with all the correct output given
//the displacement
void VertexEmit(in vec2 displacement, in float finalz)
{
  ordinate = displacement;

  vec4 shift = vec4(scale[0] * displacement, 0.0, 0.0);
  vec4 eyespace_position = gl_in[0].gl_Position + shift;
  vec4 proj_position =  ProjectionMatrix * eyespace_position;

  //Here we do the cryptic thing of inputting the finalz times by the
  //w coordinate, so that when the w divide is performed by OpenGL, we
  //will recover the w divide.
  gl_Position = vec4(proj_position.xy, finalz * proj_position.w, proj_position.w);
  EmitVertex();
}

void main()
{
  //The projected position of the center of the front of the sphere
  vec4 centerpos = ProjectionMatrix * (gl_in[0].gl_Position + vec4(0.0, 0.0, scale[0], 0.0));
  
  //We now perform the w divide so that we have the minimum
  //screen-space depth of the sphere.
  float finalz = centerpos.z / centerpos.w;

  //Standard data for each fragment
  frag_scale = scale[0];
  frag_depth = depth[0];
  model_position_frag = gl_in[0].gl_Position.xyz;
  VertexEmit(vec2(-1.0, -1.0), finalz);
  VertexEmit(vec2(-1.0, +1.0), finalz);
  VertexEmit(vec2(+1.0, -1.0), finalz);
  VertexEmit(vec2(+1.0, +1.0), finalz);
  EndPrimitive();
}
);
	}

	virtual std::string initFragmentShaderSource()
	{
	  return "#version 330\n"
	    "#ifdef GL_ARB_conservative_depth\n"
	    "#extension GL_ARB_conservative_depth : enable\n"
	    "layout (depth_greater) out float gl_FragDepth;"
	    "#endif\n"
	    STRINGIFY(
uniform mat4 ProjectionMatrix;
uniform mat4 ViewMatrix;

flat in float frag_scale;
flat in float frag_depth;
flat in vec3 model_position_frag;
smooth in vec2 ordinate;

layout (location = 0) out vec4 color_out;

void main()
{
  //The ordinate variable contains the x and y position of the
  //sphere. Use the equation of a sphere to determine the z pos
  float z = 1.0 - dot(ordinate,ordinate);

  //Discard the fragment if it lies outside the sphere
  if (z <= 0.0) discard;

  //Calculate the fragments real position on the sphere
  z = sqrt(z);
  vec4 frag_position_eye = vec4(model_position_frag + vec3(ordinate, z) * frag_scale, 1.0);
  //Calculate the fragments depth
  vec4 pos = ProjectionMatrix * frag_position_eye;
  gl_FragDepth = (pos.z / pos.w + 1.0) / 2.0;

  float depth = frag_depth;
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
