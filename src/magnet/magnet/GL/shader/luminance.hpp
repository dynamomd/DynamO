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
  L = max(10.0e-8, L);
  L_out = vec4(log(L), 1.0, 1.0, 1.0);
});
	}
      };

      class LuminanceMipMapShader: public detail::SSShader
      {
      public:
	virtual std::string initFragmentShaderSource()
	{
	  return "#version 330\n"
	    STRINGIFY(
layout (location = 0) out vec4 L_out;
smooth in vec2 screenCoord;

uniform sampler2D luminanceTex;
uniform ivec2 oldDimensions;
uniform vec2 oldInvDimensions;

vec4 data = vec4(0.0, 0.0, 0.0, 0.0);
float divider = 0.0;

vec2 oldPixelOrigin;

void operation(in ivec2 offset)
{
  //Fetch the local luminance avg (log), and maximum.
  vec2 sample = textureOffset(luminanceTex, oldPixelOrigin, offset).rg;

  //Store the value for averaging
  data.r += sample.r;
  divider += 1.0;

  //Store the maximum value
  data.g = max(sample.g, data.g);
}

void main()
{
  oldPixelOrigin = (2.0 * gl_FragCoord.xy - vec2(0.5, 0.5)) * oldInvDimensions;

  //This is the texture coordinates of the center of the lower left
  //pixel to be sampled. This is the "origin" pixel and we are going
  //to sum up the pixels above and to the right of this pixel.
  //oldPixelOrigin = screenCoord - 0.0 * oldInvDimensions;

  //First sample the standard 2x2 grid of pixels
  operation(ivec2(0,0));
  operation(ivec2(0,1));
  operation(ivec2(1,0));
  operation(ivec2(1,1));

  //Now determine if we need to add extra samples in case of
  //non-power of two textures
  bool extraXSamples = (2 * (int(gl_FragCoord.x) + 1) == oldDimensions.x - 1);
  bool extraYSamples = (2 * (int(gl_FragCoord.y) + 1) == oldDimensions.y - 1);
  
  if (extraXSamples)
    {
      operation(ivec2(2,0));
      operation(ivec2(2,1));
    }
    
  if (extraYSamples)
    {
      operation(ivec2(0,2));
      operation(ivec2(1,2));
    }

  if (extraXSamples && extraYSamples)
    operation(ivec2(2,2));

  L_out = vec4(data.r / divider, data.g, 1.0, 1.0);
});
	}
      };

    }
  }
}

#undef STRINGIFY
