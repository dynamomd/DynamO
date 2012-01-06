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
      namespace detail {
	/*! \brief A base class for shaders which downsample a texture
	  to another texture 1/4 of its resolution.
      
	  Each derived downsampling shader needs to reimplement the
	  \ref glsl_operation function.
	*/
	class DownsamplerShader: public SSShader 
	{
	  /*! \brief This function defines the two GLSL functions used
	    to combine and output the result of the samples.
	  
	    Two functions (and any global variables required) must be
	    defined, the combine and the output function.

	    The combine function is used to combine a sample to the
	    output value. This function is usually called 4 times per
	    output fragment, but may be called up to 9 times for
	    border pixels in NPOT input textures.

	    The output function is called at the end of the fragment
	    shader and must generate the value to be outputted for the
	    fragment.
	    
	    An example implementation of a averaging downsampler is
	    given below:
	    \code
vec4 sum(0.0, 0.0, 0.0, 0.0);
float counter = 0.0;

void combine(in vec4 sample)
{
  sum += sample;
  counter += 1.0;
}

vec4 output()
{
  return sum / counter;
}
	    \endcode
	  */
	  virtual std::string glsl_operation() = 0;

	  virtual std::string initFragmentShaderSource()
	  {
	    return std::string("#version 330\n") +
	      STRINGIFY(
layout (location = 0) out vec4 L_out;
smooth in vec2 screenCoord;

uniform sampler2D luminanceTex;
uniform ivec2 oldDimensions;
uniform vec2 oldInvDimensions;
) 
	      + glsl_operation() 
	      + STRINGIFY(
void main()
{
  //This is the texture coordinates of the center of the lower left
  //pixel to be sampled. This is the "origin" pixel and we are going
  //to sum up the pixels above and to the right of this pixel.
  vec2 oldPixelOrigin = (2.0 * gl_FragCoord.xy - vec2(0.5, 0.5)) * oldInvDimensions;

  //First sample the standard 2x2 grid of pixels
  combine(textureOffset(luminanceTex, oldPixelOrigin, ivec2(0,0)).rg);
  combine(textureOffset(luminanceTex, oldPixelOrigin, ivec2(0,1)).rg);
  combine(textureOffset(luminanceTex, oldPixelOrigin, ivec2(1,0)).rg);
  combine(textureOffset(luminanceTex, oldPixelOrigin, ivec2(1,1)).rg);

  //Now determine if we need to add extra samples in case of
  //non-power of two textures
  bool extraXSamples = (2 * (int(gl_FragCoord.x) + 1) == oldDimensions.x - 1);
  bool extraYSamples = (2 * (int(gl_FragCoord.y) + 1) == oldDimensions.y - 1);
  
  if (extraXSamples)
    {
      combine(textureOffset(luminanceTex, oldPixelOrigin, ivec2(2,0)).rg);
      combine(textureOffset(luminanceTex, oldPixelOrigin, ivec2(2,1)).rg);
    }
    
  if (extraYSamples)
    {
      combine(textureOffset(luminanceTex, oldPixelOrigin, ivec2(0,2)).rg);
      combine(textureOffset(luminanceTex, oldPixelOrigin, ivec2(1,2)).rg);
    }
  
  if (extraXSamples && extraYSamples)
    combine(textureOffset(luminanceTex, oldPixelOrigin, ivec2(2,2)).rg);

  L_out = output();
});
	  }
        };
      }
    }
  }
}

#undef STRINGIFY
