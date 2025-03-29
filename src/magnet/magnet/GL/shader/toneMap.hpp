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
#include <magnet/GL/shader/detail/ssshader.hpp>
#define STRINGIFY(A) #A

namespace magnet {
namespace GL {
namespace shader {
namespace detail {
/*! \brief Useful GLSL color functions for tone mapping. */
inline std::string common_glsl_color_functions() {
  return STRINGIFY(
      vec3 RGBtoXYZ(vec3 RGB) {
        const mat3 RGB2XYZ =
            mat3(0.5141364, 0.3238786, 0.16036376, 0.265068, 0.67023428,
                 0.06409157, 0.0241188, 0.1228178, 0.84442666);
        return RGB2XYZ * RGB;
      }

      vec3 XYZtoRGB(vec3 XYZ) {
        const mat3 XYZ2RGB = mat3(2.5651, -1.1665, -0.3986, -1.0217, 1.9777,
                                  0.0439, 0.0753, -0.2543, 1.1892);
        return XYZ2RGB * XYZ;
      }

      vec3 XYZtoYxy(vec3 XYZ) {
        float invXYZsum = 1.0 / dot(XYZ, vec3(1.0));
        return vec3(XYZ.g, XYZ.x * invXYZsum, XYZ.y * invXYZsum);
      }

      vec3 YxytoXYZ(vec3 Yxy) {
        return vec3(Yxy.r * Yxy.g / Yxy.b, Yxy.r,
                    Yxy.r * (1 - Yxy.g - Yxy.b) / Yxy.b);
      }

      vec3 RGBtoYxy(vec3 RGB) { return XYZtoYxy(RGBtoXYZ(RGB)); }

      vec3 YxytoRGB(vec3 Yxy) { return XYZtoRGB(YxytoXYZ(Yxy)); }

      void toneMapLuminance(inout float L, float scene_key, float inv_avg_luma,
                            float Lpwhite) {
        // Map average luminance to the middlegrey zone by scaling pixel
        // luminance
        float Lp = L * scene_key * inv_avg_luma;
        // Compress the luminance in [0,\infty) to [0,1)
        // This is Reinhard's modified mapping with controlled burnout
        L = Lp * (1.0 + Lp / (Lpwhite * Lpwhite)) / (1.0 + Lp);
      }

      void toneMapLuminance(inout float L, float scene_key, float inv_avg_luma,
                            float Lpwhite, float cutout, float O) {
        // Map average luminance to the middlegrey zone by scaling pixel
        // luminance
        float Lp = L * scene_key * inv_avg_luma;
        // Compress the luminance in [0,\infty) to [0,1)
        // This is Reinhard's modified mapping with controlled burnout
        L = max(Lp * (1.0 + Lp / (Lpwhite * Lpwhite)) - cutout, 0.0) / (O + Lp);
      }

      vec3 toneMapRGB(vec3 RGB, float scene_key, float inv_avg_luma,
                      float Lpwhite) {
        vec3 Yxy = RGBtoYxy(RGB);
        toneMapLuminance(Yxy.r, scene_key, inv_avg_luma, Lpwhite);
        return YxytoRGB(Yxy);
      }

      vec3 gammaRGBCorrection(vec3 RGB) {
        vec3 conditional = vec3(lessThan(RGB, vec3(0.00304)));
        return (conditional * 12.92 * RGB) +
               (vec3(1.0) - conditional) *
                   pow(1.055 * RGB, vec3(1.0 / 2.4) - vec3(0.055));
      });
}
} // namespace detail

/*! \brief Tone mapping shader for HDR rendering.

This is Reinhard's simplified mapping
\code
Yxy.r = Lp / (1.0 + Lp);
\endcode

This is the Filmic mapping to preserve the crispiness of the
blacks. Note: Do not use gamma correction on this function!

\code
Lp = max(0.0, Lp - 0.004);
Yxy.r = (Lp * (6.2 * Lp + 0.5))/(Lp * (6.2 * Lp + 1.7) + 0.06);
\endcode

Simplest gamma correction
\code
vec3 gamma_rgb = pow(linear_rgb, vec3(1.0/2.2));
\endcode

Better gamma correction
\code
vec3 gamma_rgb = pow(1.055 * linear_rgb, vec3(1.0/2.4) - vec3(0.055));
\endcode

 */
class ToneMapShader : public detail::SSShader {
public:
  virtual std::string initGeometryShaderSource() {
    return STRINGIFY(layout(points) in; layout(triangle_strip) out;
                     layout(max_vertices = 3) out;

                     uniform sampler2D logLuma; uniform float scene_key;
                     uniform float Lpwhite;

                     flat out float inv_log_avg_luma;
                     flat out float frag_scene_key; flat out float frag_Lpwhite;
                     smooth out vec2 screenCoord;

                     void main() {
                       // The luminance data sampled from the smallest mipmap
                       // (1x1). This holds the average logarithm of the scene
                       // luminance in the red channel and the maximum scene
                       // luminance in the green channel.
                       vec4 luma_data = texelFetch(logLuma, ivec2(0, 0), 0);

                       // Convert the average logarithm of the scene luminance
                       // into the inverse geometric mean of the luminance. We
                       // use the inverse to save doing a division in the
                       // fragment shader.

                       float Lavg = exp(luma_data.r);
                       float Lmax = luma_data.g;
                       float Lmin = luma_data.b;
                       float invlogavgluma = 1.0 / Lavg;
                       inv_log_avg_luma = invlogavgluma;

                       // Direct pass through of scene parameters
                       frag_scene_key = scene_key;
                       frag_Lpwhite = Lpwhite;

                       // This automatic scene parameter determination was taken
                       // from "Parameter Estimation for Photographic Tone
                       // Reproduction" by Erik Reinhard. frag_scene_key = 0.18
                       // *
                       // pow(4.0, (2.0 * log2(Lavg) - log2(Lmin) - log2(Lmax))
                       // / (log2(Lmax) - log2(Lmin))); frag_Lpwhite = 1.5 *
                       // pow(2.0, log2(Lmax) - log2(Lmin) - 5.0);

                       // Here we draw a fullscreen triangle and allow the GPU
                       // to scissor to the screen. This prevents the difficult
                       // interpolation of the vertex properties (e.g.,
                       // screenCoord) on the diagonal of a fullscreen quad
                       screenCoord = vec2(0.0, 0.0);
                       gl_Position = vec4(-1.0, -1.0, 0.5, 1.0);
                       EmitVertex();

                       screenCoord = vec2(2.0, 0.0);
                       gl_Position = vec4(+3.0, -1.0, 0.5, 1.0);
                       EmitVertex();

                       screenCoord = vec2(0.0, 2.0);
                       gl_Position = vec4(-1.0, +3.0, 0.5, 1.0);
                       EmitVertex();
                       EndPrimitive();
                     });
  }

  virtual std::string initFragmentShaderSource() {
    return detail::common_glsl_color_functions() +
           STRINGIFY(
               layout(location = 0) out vec4 color_out;

               uniform sampler2D color_tex; uniform sampler2D logLuma;

               uniform sampler2D bloom_tex; uniform bool bloom_enable;
               uniform float bloomCompression; uniform float bloomCutoff;
               uniform vec3 background_color;

               flat in float inv_log_avg_luma; flat in float frag_scene_key;
               flat in float frag_Lpwhite; smooth in vec2 screenCoord;

               // Taken from
               // http://www.gamedev.net/topic/407348-reinhards-tone-mapping-operator/
               void main() {
                 vec4 final_color = vec4(background_color, 0.0);

                 // Grab the scene color and tone map it
                 vec4 scene_color =
                     texelFetch(color_tex, ivec2(gl_FragCoord.xy), 0);
                 if (scene_color.a != 0.0)
                   final_color =
                       vec4(mix(final_color.rgb,
                                toneMapRGB(scene_color.rgb, frag_scene_key,
                                           inv_log_avg_luma, frag_Lpwhite),
                                scene_color.a),
                            scene_color.a);

                 // Test if bloom is enabled
                 if (bloom_enable) {
                   // Grab the blurred color, and tonemap the bloom/glare
                   vec3 bloom_RGB = texture(bloom_tex, screenCoord, 0).rgb;
                   vec3 bloom_Yxy = RGBtoYxy(bloom_RGB);
                   toneMapLuminance(bloom_Yxy.r, frag_scene_key,
                                    inv_log_avg_luma, frag_Lpwhite, bloomCutoff,
                                    1.0 / (bloomCompression + 1.0e-8));
                   bloom_RGB = YxytoRGB(bloom_Yxy);
                   final_color.rgb += max(vec3(0.0), bloom_RGB);
                 }

                 // Finally, gamma correct the image and output
                 color_out =
                     vec4(gammaRGBCorrection(final_color.rgb), final_color.a);
               });
  }
};
} // namespace shader
} // namespace GL
} // namespace magnet

#undef STRINGIFY
