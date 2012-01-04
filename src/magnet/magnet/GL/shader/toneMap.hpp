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
	virtual std::string initVertexShaderSource()
	{ return "#version 330\n"
	    STRINGIFY(
layout (location = 0) in vec4 vPosition;
uniform sampler2D logLuma;
flat out float inv_log_avg_luma;
uniform sampler2D color_tex;

void main() 
{
  const vec2 madd=vec2(0.5, 0.5);
  gl_Position = vec4(vPosition.xy, 0.0, 1.0);
  inv_log_avg_luma = 1.0 / exp(textureLod(logLuma, vec2(0.5, 0.5), 100.0).r);
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
uniform float invGamma;
uniform float exposure;

flat in float inv_log_avg_luma;

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
  
  //Compress the luminance in [0,\infty)
  Yxy.r = Lp / (1.0 + Lp);

  //Convert from Yyx to XYZ
  XYZ = vec3(Yxy.r * Yxy.g / Yxy.b, Yxy.r, Yxy.r * (1 - Yxy.g - Yxy.b) / Yxy.b);

  //Convert from XYZ to RGB
  const mat3 XYZ2RGB = mat3(2.5651,-1.1665,-0.3986,
			    -1.0217, 1.9777, 0.0439, 
			    0.0753, -0.2543, 1.1892);
  color_out = vec4(XYZ2RGB * XYZ, 1.0);
});
	}
      };
    }
  }
}

#undef STRINGIFY
