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
      /*! \brief A Screen Space Ambient Occlusion Shader.
       *
       * This effect is used to fake the ambient light effects caused
       * by indirect lighting. There are no good tutorials on this
       * method, but here are some links to some more information.
       *
       * http://www.gamerendering.com/category/lighting/ssao-lighting/
       */
      class SSAOShader : public detail::SSShader
      {
      public:
	virtual std::string initFragmentShaderSource()
	{
	  return STRINGIFY(
uniform sampler2DMS NormalsTex;
uniform sampler2DMS EyePosTex;
uniform sampler2D rnm;
uniform float radius;
uniform float totStrength;
uniform float depthDropoff;
uniform float offset;

uniform mat4 ProjectionMatrix;
uniform mat4 ViewMatrix;

const float invSamples = 1.0 / 10.0;

smooth in vec2 screenCoord;
layout (location = 0) out vec4 color_out;

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
  vec3 fres = normalize(2.0 * texture(rnm, screenCoord * offset).xyz - vec3(1.0));
    
  ivec2 pixelcoord = ivec2(textureSize(EyePosTex) * screenCoord);
  vec3 currentPixelPos = texelFetch(EyePosTex, pixelcoord, 0).xyz; 
  vec3 currentPixelNorm = normalize(2.0 * texelFetch(NormalsTex, pixelcoord, 0).xyz - 1.0);

  float occlusion = 0.0;
  for(int i = 0; i < 10; ++i)
    {
      //Get a vector (randomized inside of a sphere with radius 1.0)
      //from a texture and reflect it through the random vector.
      vec3 ray = radius * reflect(pSphere[i], fres);
      
      //Calculate the eye space location of the sample, we make sure
      //the sample is in the correct hemisphere around the normal
      //using the sign function.
      vec3 sampleEyeSpaceLocation =  currentPixelPos + sign(dot(ray, currentPixelNorm)) * ray;

      //Convert it to clip space
      vec4 sampleLocation = ProjectionMatrix * (ViewMatrix * vec4(sampleEyeSpaceLocation, 1.0));
      
      //Convert to normalised device coordinates
      sampleLocation.xyz *= (1.0 / sampleLocation.w);
      
      //Change to normalised screen coordinates
      sampleLocation.xyz = 0.5 * sampleLocation.xyz + vec3(1.0);

      //Now fetch the scenes normal and position at the screen
      //position of the sample. 
      ivec2 samplecoord = ivec2(textureSize(EyePosTex) * sampleLocation.xy);
      vec3 occluderLocation = texelFetch(EyePosTex, samplecoord, 0).xyz;
      vec3 occluderNormal = normalize(texelFetch(NormalsTex, samplecoord, 0).xyz * 2.0 - vec3(1.0));

      //The z coordinate is increasing as objects get farther away, so
      //if the zdiff is negative the sample is behind the occluder.
      float zdiff =  occluderLocation.z - sampleEyeSpaceLocation.z;

      //The curvature of the local surroundings is also important
      float occluderDot = dot(occluderNormal, currentPixelNorm);
      
      occlusion += step(0.0, -zdiff);
      //occlusion += max(0.0, 0.9 - occluderDot) * step(0.0, zdiff) * step(0.0, depthDropoff - zdiff);
    }

  float val = clamp(1.0 - totStrength * occlusion * invSamples, 0.0, 1.0);

  color_out = vec4(vec3(val), 1.0);
});
	}
      };
    }
  }
}

#undef STRINGIFY
