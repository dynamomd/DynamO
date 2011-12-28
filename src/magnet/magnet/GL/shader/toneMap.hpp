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
smooth out vec2 screenCoord;
uniform sampler2D logLuma;
flat out float inv_avg_luma;
uniform sampler2D color_tex;

void main() 
{
  const vec2 madd=vec2(0.5, 0.5);
  screenCoord = vPosition.xy * madd + madd;
  gl_Position = vec4(vPosition.xy, 0.0, 1.0);
  //  inv_avg_luma = 1.0 / exp(textureLod(logLuma, vec2(0.5, 0.5), 100.0).r);
  inv_avg_luma = 1.0 / dot(textureLod(color_tex, vec2(0.5, 0.5), 100.0).rgb, 
			   vec3(0.2126, 0.7152, 0.0722));
}); 
	}

	virtual std::string initFragmentShaderSource()
	{
	  return "#version 330\n"
	    STRINGIFY(
//Normalized position on the screen
smooth in vec2 screenCoord;
layout (location = 0) out vec4 color_out;

uniform sampler2D color_tex;
uniform sampler2D logLuma;
uniform float invGamma;
uniform float exposure;

flat in float inv_avg_luma;

void main()
{
  vec3 color_in = texture(color_tex, screenCoord).rgb;
  float old_luminance = dot(color_in, vec3(0.2126, 0.7152, 0.0722));
  float relative_luma = exposure * inv_avg_luma;
  //color_out = vec4(pow(relative_luma  * color_in, vec3(invGamma)), 1.0);
  //float val = exp(textureLod(logLuma, screenCoord, 0.0).r);
  float val = exp(textureLod(logLuma, screenCoord, 300.0).r);
  color_out = vec4(val,val,val, 1.0);
});
	}
      };
    }
  }
}

#undef STRINGIFY
