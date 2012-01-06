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
smooth out vec2 screenCoord;

void main()
{
  //Here we draw a fullscreen triangle and allow the GPU to scissor to
  //the screen. This prevents the difficult interpolation of the
  //vertex property (screenCoord) on the diagonal of a fullscreen
  //quad. This is a ridiculous optimisation I know.

  vec2 luma_data = textureLod(logLuma, vec2(0.5, 0.5), 100.0).rg;
  float invlogavgluma = 1.0 / exp(luma_data.r);
  inv_log_avg_luma = invlogavgluma;

  float tmp = burnout * luma_data.g * exposure * invlogavgluma;
  scaled_burnout_luma = tmp * tmp;

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
//Normalized position on the screen
layout (location = 0) out vec4 color_out;

uniform sampler2D color_tex;
uniform sampler2D logLuma;
uniform float exposure;

flat in float inv_log_avg_luma;
flat in float scaled_burnout_luma;

//Taken from http://www.gamedev.net/topic/407348-reinhards-tone-mapping-operator/
void main()
{
  vec3 color_in = texelFetch(color_tex, ivec2(gl_FragCoord.xy), 0).rgb;

  //Convert from RGB to XYZ
  const mat3 RGB2XYZ = mat3(0.5141364, 0.3238786,  0.16036376,
			    0.265068,  0.67023428, 0.06409157,
			    0.0241188, 0.1228178,  0.84442666);
  vec3 XYZ = RGB2XYZ * color_in;
  
  //Convert from XYZ to Yxy
  float invXYZsum = 1.0 / dot(XYZ, vec3(1.0, 1.0, 1.0));
  vec3 Yxy = vec3(XYZ.g, XYZ.r * invXYZsum, XYZ.g * invXYZsum);

  // (Lp) Map average luminance to the middlegrey zone by scaling pixel luminance
  float Lp = Yxy.r * exposure * inv_log_avg_luma;
  
  //Compress the luminance in [0,\infty) to [0,1)

  //This is Reinhard's simplified mapping
  //Yxy.r = Lp / (1.0 + Lp);

  //This is Reinhard's modified mapping with controlled burnout
  Yxy.r = Lp * (1.0 + Lp / scaled_burnout_luma)  / (1.0 + Lp);

  //This is the Filmic mapping to preserve the crispiness of the
  //blacks. Note: Do not use gamma correction on this function!
  //Lp = max(0.0, Lp - 0.004); 
  //Yxy.r = (Lp * (6.2 * Lp + 0.5))/(Lp * (6.2 * Lp + 1.7) + 0.06);

  //Convert from Yyx to XYZ
  XYZ = vec3(Yxy.r * Yxy.g / Yxy.b, Yxy.r, Yxy.r * (1 - Yxy.g - Yxy.b) / Yxy.b);

  //Convert from XYZ to RGB
  const mat3 XYZ2RGB = mat3(2.5651,-1.1665,-0.3986,
			    -1.0217, 1.9777, 0.0439, 
			    0.0753, -0.2543, 1.1892);

  //Additional gamma correction
  vec3 linear_rgb = XYZ2RGB * XYZ;

  //Simplest gamma correction
  //vec3 gamma_rgb = pow(linear_rgb, vec3(1.0/2.2));

  //Better gamma correction
  //vec3 gamma_rgb = pow(1.055 * linear_rgb, vec3(1.0/2.4) - vec3(0.055));
  
  //Best gamma correction
  vec3 conditional = vec3(lessThan(linear_rgb, vec3(0.00304)));
  vec3 gamma_rgb = conditional * 12.92 * linear_rgb;
  gamma_rgb += (vec3(1.0) - conditional)
    * pow(1.055 * linear_rgb, vec3(1.0/2.4) - vec3(0.055));
  
  color_out = vec4(gamma_rgb, 1.0);
});
	}
      };
    }
  }
}

#undef STRINGIFY
