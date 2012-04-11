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
      /*! \brief Implements a Bilateral 5x5 Gaussian Blur Shader.
       *
       * A Bilateral blur is one that takes depth information into
       * account and will not blur across sharp changes in the depth.
       *
       * This is useful when trying to blur the surface of an object,
       * but to avoid blurring its edge.
       */
      class BilateralBlur : public detail::SSShader
      {
      public:
	virtual std::string initFragmentShaderSource()
	{ 
	  return STRINGIFY(
uniform sampler2D ImageTex; //input
uniform sampler2DMS EyePosTex;
uniform float totStrength;
uniform float nearDist;
uniform float farDist;
uniform int radius;

layout (location = 0) out vec4 color_out;

const float weight[5] = float[5](0.05496597,0.24581,0.4076311347,0.24581,0.05496597);

float sampleWeight(int i, int j) { return weight[i] * weight[j]; }

void main(void)
{
  float currentPixelDepth = texelFetch(EyePosTex, ivec2(gl_FragCoord.xy),0).z;
  
  vec3 accum = vec3(0, 0, 0);
  float totalWeight = 0.0;
  
  for (int x = 0; x < 5; ++x)
    for (int y = 0; y < 5; ++y)
      {
	ivec2 sample_coords = ivec2(gl_FragCoord.xy) + ivec2(x - 2, y - 2) * radius;
	float sampleDepth = texelFetch(EyePosTex, sample_coords, 0).z;
	
	float Zdifference = abs(currentPixelDepth - sampleDepth);
	float sampleweight = (1.0 - step(totStrength, Zdifference)) * sampleWeight(x,y);
	accum += sampleweight * texelFetch(ImageTex, sample_coords, 0).rgb;
	totalWeight += sampleweight;
      }
  
  color_out = vec4(accum / totalWeight, 1);
});
	}
      };
    }
  }
}

#undef STRINGIFY
