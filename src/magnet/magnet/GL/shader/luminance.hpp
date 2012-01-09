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
#include <magnet/GL/shader/downsampler.hpp>
#define STRINGIFY(A) #A

namespace magnet {
  namespace GL {
    namespace shader {
      /*! \brief This calculates the logarthmic luminance values for
          the pixels in a scene
       */
      class LuminanceShader: public detail::SSShader
      {
      public:
	virtual std::string initFragmentShaderSource()
	{
	  return "#version 330\n"
	    STRINGIFY(
//Normalized position on the screen
smooth in vec2 screenCoord;
layout (location = 0) out vec4 L_out;

//The HDR color buffer
uniform sampler2D colorTex;

void main()
{
  vec4 color = texture(colorTex, screenCoord).rgba;
  float L = dot(color.rgb, vec3(0.265068,  0.67023428, 0.06409157));
  //Prevent negative logarithms
  L_out = vec4(log(max(10.0e-8, L)), L, 1.0, 1.0);
});
	}
      };

      class LuminanceMipMapShader: public DownsamplerShader
      {
      public:
	virtual std::string glsl_operation()
	{
	  return STRINGIFY(
vec2 data = vec2(0.0);
float divider = 0.0;

void combine(in vec4 sample)
{
  //Store the value for averaging
  data.r += sample.r;
  divider += 1.0;

  //Store the maximum value
  data.g = max(sample.g, data.g);
}

vec4 output_frag()
{
  return vec4(data.r / divider, data.g, 1.0, 1.0);
}
);
	}
      };

    }
  }
}

#undef STRINGIFY
