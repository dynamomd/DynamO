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
      class SSAO : public detail::Shader
      {
      public:
	void invoke()
	{
	  attach();
	  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	  drawScreenQuad();
	  //Restore the fixed pipeline
	  glUseProgramObjectARB(0);
	}

	virtual std::string initVertexShaderSource()
	{
	  return STRINGIFY(void main(void) { gl_Position = ftransform(); gl_TexCoord[0] = gl_MultiTexCoord0; });
	}

	virtual std::string initFragmentShaderSource()
	{
	  return STRINGIFY(
uniform sampler2D u_Texture1;
uniform sampler2D u_Texture2;
uniform sampler2D rnm;
uniform float radius;
uniform float totStrength;
uniform float depthDropoff;
uniform float offset;
uniform float nearDist;
uniform float farDist;

const float invSamples = 1.0 / 10.0;

float LinearizeDepth(float zoverw)
{
  return(2.0 * nearDist) / (farDist + nearDist - zoverw * (farDist - nearDist));
}

void main(void)
{
  // these are the random vectors inside a unit sphere
  vec3 pSphere[10] = vec3[10](vec3(-0.010735935, 0.01647018, 0.0062425877),
			      vec3(-0.06533369, 0.3647007, -0.13746321),
			      vec3(-0.6539235, -0.016726388, -0.53000957),
			      vec3(0.40958285, 0.0052428036, -0.5591124),
			      vec3(-0.1465366, 0.09899267, 0.15571679),
			      vec3(-0.44122112, -0.5458797, 0.04912532),
			      vec3(0.03755566, -0.10961345, -0.33040273),
			      vec3(0.019100213, 0.29652783, 0.066237666),
			      vec3(0.8765323, 0.011236004, 0.28265962),
			      vec3(0.29264435, -0.40794238, 0.15964167));

  // grab a normal for reflecting the sample rays later on
  vec3 fres = normalize(2.0 * texture2D(rnm, gl_TexCoord[0].st * offset).xyz - vec3(1.0));
    
  float currentPixelDepth = LinearizeDepth(texture2D(u_Texture2, gl_TexCoord[0].st).r);
  
  // current fragment coords in screen space
  vec3 ep = vec3(gl_TexCoord[0].st, currentPixelDepth);
  // get the normal of current fragment
  vec3 norm = normalize(2.0 * texture2D(u_Texture1, gl_TexCoord[0].st).xyz - 1.0);
  
  float bl = 0.0;
  float radD = radius / currentPixelDepth;
  
  //vec3 ray, se, occNorm;
  float occluderDepth;
  for(int i = 0; i < 10; ++i)
    {
      // get a vector (randomized inside of a sphere with radius 1.0) from a texture and reflect it
      vec3 ray = radD * reflect(pSphere[i],fres);
      
      // get the depth of the occluder fragment
      vec2 occluderLoc = ep.xy + sign(dot(ray,norm) ) * ray.xy;
      vec4 occluderFragment = texture2D(u_Texture1, occluderLoc);

      float occluderDepth = LinearizeDepth(texture2D(u_Texture2, occluderLoc).r);

      // if d (depth difference) is negative = occluder is behind current fragment
      float d = currentPixelDepth - occluderDepth;
      
      vec3 occluderNorm = normalize(occluderFragment.xyz * 2.0 - 1.0);
      float occluderDot = dot(occluderNorm, norm);
      
      bl += max(0, 0.9 - occluderDot) * step(0, d) * step(0, depthDropoff - d);
    }

  float val = clamp(1.0 - totStrength * bl * invSamples, 0, 1);

  gl_FragColor = vec4(val, val, val, 1.0);
});
	}
      };
    }
  }
}

#undef STRINGIFY
