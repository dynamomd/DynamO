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

    A discussion of this technique is given in the excellent
    online GL book by Jason L. McKesson at \url
    http://www.arcsynthesis.org/gltut/index.html in the chapter
    on lies and IMPOSTORS.
 */
class SphereShader : public detail::Shader {
public:
  SphereShader() { defines("unshaded") = "false"; }

  virtual std::string initVertexShaderSource() {
    return STRINGIFY(uniform mat4 ViewMatrix; uniform float global_scale;

                     layout(location = 0) in vec4 vPosition;
                     layout(location = 1) in vec4 vColor;
                     layout(location = 4) in vec4 iOrientation;
                     layout(location = 5) in vec4 iScale;

                     out vec4 color; out float radius;

                     void main() {
                       color = vColor;
                       radius = (iScale.x + float(iScale.x == 0.0)) *
                                global_scale * 0.5;
                       gl_Position = ViewMatrix * vec4(vPosition.xyz, 1.0);
                     });
  }

  virtual std::string initGeometryShaderSource() {
    return STRINGIFY(
        uniform mat4 ProjectionMatrix;

        layout(points) in; layout(triangle_strip) out;
        layout(max_vertices = 4) out;

        in vec4 color[]; in float radius[];

        flat out vec4 vert_color; flat out vec3 frag_center;
        flat out float frag_radius; smooth out vec3 frag_pos;

        // Function to emit a bilboard vertex with all the correct output given
        // the displacement
        void VertexEmit(in vec2 displacement) {
          // The billboards need to be slightly larger to accommodate
          // perspective warping.
          const float overdraw = 1.1;
          displacement *= overdraw;
          frag_radius = radius[0];
          vert_color = color[0];
          frag_center = gl_in[0].gl_Position.xyz;
          vec3 position = gl_in[0].gl_Position.xyz +
                          vec3(radius[0] * displacement, radius[0]);
          frag_pos = position;
          gl_Position =
              ProjectionMatrix * vec4(position, gl_in[0].gl_Position.w);
          EmitVertex();
        }

        void main() {
          // Standard data for each fragment
          VertexEmit(vec2(-1.0, -1.0));
          VertexEmit(vec2(-1.0, +1.0));
          VertexEmit(vec2(+1.0, -1.0));
          VertexEmit(vec2(+1.0, +1.0));
          EndPrimitive();
        });
  }

  virtual std::string initFragmentShaderSource() {
    return "#ifdef GL_ARB_conservative_depth\n"
           "#extension GL_ARB_conservative_depth : enable\n"
           "layout (depth_greater) out float gl_FragDepth;"
           "#endif\n" STRINGIFY(
               uniform mat4 ProjectionMatrix;

               flat in vec4 vert_color; flat in vec3 frag_center;
               flat in float frag_radius; smooth in vec3 frag_pos;

               layout(location = 0) out vec4 color_out;
               layout(location = 1) out vec4 normal_out;
               layout(location = 2) out vec4 position_out;

               void main() {
                 normal_out = vec4(0.0);
) "\n#ifdef DRAWBILLBOARD\n" STRINGIFY(
  position_out = vec4(frag_pos, 1.0);
  vec4 screen_pos = ProjectionMatrix * vec4(frag_pos, 1.0);
  gl_FragDepth = (screen_pos.z / screen_pos.w + 1.0) / 2.0;
) "\n#else\n" STRINGIFY(
  vec3 rij = -frag_center;
  vec3 vij = frag_pos;
  gl_FragDepth = gl_FragCoord.z;
  
  float A = dot(vij, vij);
  float B = dot(rij, vij);
  float C = dot(rij, rij) - frag_radius * frag_radius;
  float argument = B * B - A * C;
  if (argument < 0.0) discard;
  float t = - C / (B - sqrt(argument));
  vec3 hit = t * vij;
  position_out = vec4(hit, 1.0);
  
  if (!unshaded)
    normal_out = vec4(normalize(hit - frag_center),1.0);

  vec4 screen_pos = ProjectionMatrix * vec4(hit, 1.0);
  gl_FragDepth = (screen_pos.z / screen_pos.w + 1.0) / 2.0;
) "\n#endif\n" STRINGIFY(
  color_out = vert_color;
               });
  }
};

/*! \brief A variant of the SphereShader used for variance
    shadow mapping.
 */
class SphereVSMShader : public SphereShader {
  virtual std::string initFragmentShaderSource() {
    return STRINGIFY(uniform mat4 ProjectionMatrix;

                     flat in vec3 frag_center; flat in float frag_radius;
                     smooth in vec3 frag_pos;

                     layout(location = 0) out vec4 moments_out;

                     void main() {
) "\n#ifdef DRAWBILLBOARD\n" STRINGIFY(
  vec3 world_pos = frag_pos;
) "\n#else\n" STRINGIFY(
  vec3 rij = -frag_center;
  vec3 vij = frag_pos;
  
  float A = dot(vij, vij);
  float B = dot(rij, vij);
  float C = dot(rij, rij) - frag_radius * frag_radius;
  float argument = B * B - A * C;
  if (argument < 0.0) discard;
  float t = - C / (B - sqrt(argument));
  vec3 world_pos = t * vij;
) "\n#endif\n" STRINGIFY(
  vec4 screen_pos = ProjectionMatrix * vec4(world_pos, 1.0);
  gl_FragDepth = (screen_pos.z / screen_pos.w + 1.0) / 2.0;

  float moment1 = length(world_pos);
  float moment2 = moment1 * moment1;
  float dx = dFdx(moment1);
  float dy = dFdy(moment1);
  moment2 += 0.25 * (dx * dx + dy * dy);
  moments_out = vec4(moment1, moment2, 0.0, 1.0);
                     });
  }
};
} // namespace shader
} // namespace GL
} // namespace magnet

#undef STRINGIFY
