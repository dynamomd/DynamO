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
#include <magnet/GL/shader/detail/ssshader.hpp>
#define STRINGIFY(A) #A

namespace magnet {
  namespace GL {
    namespace shader {
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
      class ToneMapShader: public detail::SSShader
      {
      public:
	virtual std::string initGeometryShaderSource()
	{ return "#version 330\n"
	    STRINGIFY(
layout(points) in;
layout(triangle_strip) out;
layout(max_vertices = 3) out;

uniform sampler2D logLuma;
uniform float exposure;
uniform float burnout;

flat out float inv_log_avg_luma;
flat out float scaled_burnout_luma;
flat out float burnout_luma;
smooth out vec2 screenCoord;

void main()
{
  //Here we draw a fullscreen triangle and allow the GPU to scissor to
  //the screen. This prevents the difficult interpolation of the
  //vertex property (screenCoord) on the diagonal of a fullscreen
  //quad. This is a ridiculous optimisation I know.

  vec2 luma_data = textureLod(logLuma, vec2(0.5, 0.5), 100.0).rg;
  float log_avg_luma = exp(luma_data.r);
  float invlogavgluma = 1.0 / log_avg_luma;
  inv_log_avg_luma = invlogavgluma;

  burnout_luma = burnout * luma_data.g;

  float tmp = burnout * luma_data.g * exposure * invlogavgluma;
  scaled_burnout_luma = tmp;

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

	virtual std::string initFragmentShaderSource()
	{
	  return "#version 330\n"
	    STRINGIFY(
layout (location = 0) out vec4 color_out;

uniform sampler2D color_tex;
uniform sampler2D logLuma;
uniform float exposure;

uniform sampler2D bloom_tex;
uniform float bloomStrength;
uniform float bloomCutoff;

flat in float inv_log_avg_luma;
flat in float scaled_burnout_luma;
flat in float burnout_luma;
smooth in vec2 screenCoord;


vec3 RGBtoXYZ(vec3 input)
{
  const mat3 RGB2XYZ = mat3(0.5141364, 0.3238786,  0.16036376,
			    0.265068,  0.67023428, 0.06409157,
			    0.0241188, 0.1228178,  0.84442666);
  return RGB2XYZ * input;
}

vec3 XYZtoRGB(vec3 input)
{
  const mat3 XYZ2RGB = mat3(2.5651,-1.1665,-0.3986,
			    -1.0217, 1.9777, 0.0439, 
			    0.0753, -0.2543, 1.1892);
  return XYZ2RGB * input;
}

vec3 XYZtoYxy(vec3 input)
{
  float invXYZsum = 1.0 / dot(input, vec3(1.0));
  return vec3(input.g, input.r * invXYZsum, input.g * invXYZsum);
}

vec3 YxytoXYZ(vec3 input)
{
  return vec3(input.r * input.g / input.b, 
	      input.r,
	      input.r * (1 - input.g - input.b) / input.b);
}

vec3 RGBtoYxy(vec3 input)
{ return XYZtoYxy(RGBtoXYZ(input)); }

vec3 YxytoRGB(vec3 input)
{ return XYZtoRGB(YxytoXYZ(input)); }

vec3 toneMapYxy(vec3 input, float exposure, 
		float inv_avg_luma, 
		float scaled_burnout_luma)
{
  //Map average luminance to the middlegrey zone by scaling pixel luminance
  float Lp = input.r * exposure * inv_avg_luma;

  scaled_burnout_luma *= scaled_burnout_luma;
  
  //Compress the luminance in [0,\infty) to [0,1)
  //This is Reinhard's modified mapping with controlled burnout
  input.r = Lp * (1.0 + Lp / scaled_burnout_luma) / (1.0 + Lp);

  return input;
}

vec3 gammaRGBCorrection(vec3 input)
{
  vec3 conditional = vec3(lessThan(input, vec3(0.00304)));
  return (conditional * 12.92 * input)
    + (vec3(1.0) - conditional)
    * pow(1.055 * input, vec3(1.0/2.4) - vec3(0.055));
}

//Taken from http://www.gamedev.net/topic/407348-reinhards-tone-mapping-operator/
void main()
{
  vec3 bloom_RGB = texture(bloom_tex, screenCoord, 0).rgb;
  vec3 bloom_Yxy = RGBtoYxy(bloom_RGB);
  bloom_Yxy.r = bloomStrength * max(bloom_Yxy.r - bloomCutoff * burnout_luma, 0.0);
  

  vec3 color_RGB = texelFetch(color_tex, ivec2(gl_FragCoord.xy), 0).rgb;

  color_RGB += YxytoRGB(bloom_Yxy);

  vec3 scaled_RGB = YxytoRGB(toneMapYxy(RGBtoYxy(color_RGB), 
					exposure,
					inv_log_avg_luma, 
					scaled_burnout_luma));

  vec3 gamma_RGB = gammaRGBCorrection(scaled_RGB);
  color_out = vec4(gamma_RGB, 1.0);
});
	}
      };
    }
  }
}

#undef STRINGIFY
