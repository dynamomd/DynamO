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
      /*! \brief A base class for shaders which downsample a texture
	to another texture 1/4 of its resolution.
      
	Each derived downsampling shader needs to reimplement the
	\ref glsl_operation function. This base class is an
	implementation of the trivial averaging shader.
      */
      class DownsamplerShader: public detail::SSShader
      {
	/*! \brief This function defines the two GLSL functions used
	  to combine and output the result of the samples.
	  
	  Two functions (and any global variables required) must be
	  defined, the combine and the output_frag function.

	  The combine function is used to combine a sample to the
	  output value. This function is usually called 4 times per
	  output fragment, but may be called up to 9 times for border
	  pixels in NPOT input textures.

	  The output_frag function is called at the end of the
	  fragment shader and must generate the value to be outputted
	  for the fragment.
	*/
	virtual std::string glsl_operation()
	{
	  return STRINGIFY(
vec4 sum = vec4(0.0);
float counter = 0.0;

void combine(in vec4 sample)
{
  sum += sample;
  counter += 1.0;
}

vec4 output_frag()
{
  return sum / counter;
}
);
	}

	virtual std::string initFragmentShaderSource()
	{
	  return STRINGIFY(
layout (location = 0) out vec4 L_out;
uniform sampler2D inputTex;
uniform ivec2 oldSize;
uniform int downscale = 2;
) 
	    + glsl_operation() 
	    + STRINGIFY(
void main()
{
  //This is the texture coordinates of the center of the lower left
  //pixel to be sampled. This is the "origin" pixel and we are going
  //to sum up the pixels above and to the right of this pixel.
  ivec2 oldPixelOrigin = downscale * ivec2(gl_FragCoord.xy);

  int step = downscale / 2;
  //First sample the standard 2x2 grid of pixels
  combine(texelFetch(inputTex, oldPixelOrigin + step * ivec2(0,0), 0));
  combine(texelFetch(inputTex, oldPixelOrigin + step * ivec2(0,1), 0));
  combine(texelFetch(inputTex, oldPixelOrigin + step * ivec2(1,0), 0));
  combine(texelFetch(inputTex, oldPixelOrigin + step * ivec2(1,1), 0));

  //Now determine if we need to add extra samples in case of
  //non-power of two textures
  bool extraXSamples = oldPixelOrigin.x + downscale == oldSize.x - 1;
  bool extraYSamples = oldPixelOrigin.y + downscale == oldSize.y - 1;
  
  if (extraXSamples)
    {
      combine(texelFetch(inputTex, oldPixelOrigin + step * ivec2(2,0), 0));
      combine(texelFetch(inputTex, oldPixelOrigin + step * ivec2(2,1), 0));
    }
    
  if (extraYSamples)
    {
      combine(texelFetch(inputTex, oldPixelOrigin + step * ivec2(0,2), 0));
      combine(texelFetch(inputTex, oldPixelOrigin + step * ivec2(1,2), 0));
    }
  
  if (extraXSamples && extraYSamples)
    combine(texelFetch(inputTex, oldPixelOrigin + step * ivec2(2,2), 0));

  L_out = output_frag();
});
	}
      };
    }
  }
}

#undef STRINGIFY
