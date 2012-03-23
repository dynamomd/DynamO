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
  vec4 color = texture(colorTex, screenCoord);
  float L = dot(color.rgb, vec3(0.265068,  0.67023428, 0.06409157));
  //Prevent logarithms of zero, store the log(L), max L, min L, weight/alpha
  L_out = vec4(log(max(1.0e-5, L)), L, L, color.a/10000.0);
  //The weight is divided by 10000.0 to use most of the range of the
  //exponent in the half-precision floating point format. (there may
  //be more than 65504 fragments in an image, but this is the max
  //16-bit floating point value. The smallest floating point value is 2^{-14}
});
	}
      };

      class LuminanceMipMapShader: public DownsamplerShader
      {
      public:
	virtual std::string glsl_operation()
	{
	  return STRINGIFY(
vec4 data = vec4(0.0);

void combine(in vec4 sample)
{
  if (sample.a != 0.0)
    {
      //If this is the first sample, just copy the min max values.
      if (data.a == 0.0)
	data.gb = sample.gb;

      //Store the value for averaging, weighted by the rendered
      //fragment count
      data.r += sample.r * sample.a;
      //Store the maximum value
      data.g = max(sample.g, data.g);
      //Store the maximum value
      data.b = min(sample.b, data.b);
      //Add on the fragment count of this sample
      data.a += sample.a;
    }
}

vec4 output_frag() { 
  if (data.a != 0.0) 
    data.r /= data.a;
  return data; 
}
);
	}
      };

    }
  }
}

#undef STRINGIFY
