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
#include <magnet/GL/shader/detail/shader.hpp>
#define STRINGIFY(A) #A

namespace magnet {
  namespace GL {
    namespace shader {
      class BilateralBlur : public detail::Shader
      {
      public:
	void invoke()
	{
	  attach();
	  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	  drawScreenQuad();
	  glUseProgramObjectARB(0);	
	}
	
	virtual std::string initFragmentShaderSource()
	{ 
	  return STRINGIFY(
uniform sampler2D u_Texture0; //input
uniform sampler2D u_Texture2; //Depth buffer
uniform vec2 scale;
uniform float totStrength;
uniform float nearDist;
uniform float farDist;
const float invSamples = 1.0 / 10.0;

const float weight[5] = float[5](0.05496597,0.24581,0.4076311347,0.24581,0.05496597);

float sampleWeight(int i, int j) { return weight[i] * weight[j]; }

float LinearizeDepth(float zoverw)
{
  return(2.0 * nearDist) / (farDist + nearDist - zoverw * (farDist - nearDist));
}

void main(void)
{
  float currentPixelDepth = LinearizeDepth(texture2D(u_Texture2, gl_TexCoord[0].st).r);
  
  vec3 accum = vec3(0, 0, 0);
  float totalWeight = 0.0;
  
  for (int x = 0; x < 5; ++x)
    for (int y = 0; y < 5; ++y)
      {
	vec2 sampleLoc = gl_TexCoord[0].st + vec2((x - 2) * scale.x, (y - 2) * scale.y);
	float sampleDepth = LinearizeDepth(texture2D(u_Texture2, sampleLoc).r);
	
	float Zdifference = abs(currentPixelDepth - sampleDepth);
	float sampleweight = (1.0 - step(totStrength, Zdifference)) * sampleWeight(x,y);
	accum += sampleweight * texture2D(u_Texture0, sampleLoc).rgb;
	totalWeight += sampleweight;
      }
  
  gl_FragColor = vec4(accum / totalWeight, 1);
});
	}
      };
    }
  }
}

#undef STRINGIFY
