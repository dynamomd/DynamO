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
      /*! \brief Deffered lighting calculation shader.

	This class performs the lighting calculations for the current
	scene.
       */
      class AmbientLightShader: public detail::SSShader
      {
      public:
	virtual std::string initFragmentShaderSource()
	{
	  return STRINGIFY(
layout (location = 0) out vec4 color_out;

//Standard G-buffer data
uniform sampler2DMS colorTex;
uniform int samples;
uniform float ambientLight;

void main()
{
  //Now calculate the color from the samples
  vec4 color_sum = vec4(0.0);
  
  for (int sample_id = 0; sample_id < samples; sample_id++)
    {
      vec4 color = texelFetch(colorTex, ivec2(gl_FragCoord.xy), sample_id).rgba;

      //If alpha is zero, this is an empty pixel, and should not
      //contribute to the tone mapping
      if (color.a != 0)
	{
	  color_sum.rgb += ambientLight * color.rgb;
	  color_sum.a += 1.0;
	}
    }
 
  //We write out the HDR color here, along with the occupancy
  //(fraction of drawn pixels) in the alpha channel.
  color_out = color_sum / float(samples);
});
	}
      };
    }
  }
}

#undef STRINGIFY
