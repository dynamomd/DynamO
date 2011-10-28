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
      /*! \brief Deffered lighting calculation shader.

	This class performs the lighting calculations for the current
	scene.
       */
      class DeferredLightingShader: public detail::SSShader
      {
      public:
	virtual std::string initFragmentShaderSource()
	{
	  return "#version 330\n"
	    STRINGIFY(
//Normalized position on the screen
smooth in vec2 screenCoord;
layout (location = 0) out vec4 color_out;

//Standard G-buffer data
uniform sampler2D colorTex;
uniform sampler2D normalTex;
uniform sampler2D positionTex;
uniform vec3 lightPosition;

///////////////Shadow mapping functions and variables
uniform mat4 ShadowMatrix;

vec4 ShadowCoord;
uniform sampler2D ShadowMap;
uniform int ShadowMapping;
uniform float ShadowIntensity;

float linstep(float min, float max, float v)  
{  
  return clamp((v - min) / (max - min), 0.0, 1.0);  
}  

float ReduceLightBleeding(float p_max, float Amount)  
{  
  // Remove the [0, Amount] tail and linearly rescale (Amount, 1].  
  return linstep(Amount, 1.0, p_max);  
}  

float chebyshevUpperBound(float distance)
{
  vec2 moments = texture2D(ShadowMap,ShadowCoord.xy).rg;
	
  // We use chebyshev's upperBound to check How likely this pixel is
  // to be lit (p_max)
  float variance = moments.y - (moments.x * moments.x);
  variance = max(variance, 0.0000001);

  float d = distance - moments.x;
  float p_max = variance / (variance + d * d);

  //We linearly remap the probability so that a certain range is
  //always completely in shadow
  p_max = ReduceLightBleeding(p_max, 0.2);

  float p = float(distance <= moments.x);
  return max(p, p_max);
}

void main()
{
  //Eye space position of the vertex
  vec3 position = texture2D(positionTex, screenCoord).xyz;
  
  //Eye space normal of the vertex
  vec3 normal = texture2D(normalTex, screenCoord).rgb;
  normal = normalize(normal);

  //light position relative to the pixel location
  vec3 lightVector = lightPosition - position;
  float lightDistance = length(lightVector);
  vec3 lightDirection = lightVector * (1.0 / lightDistance);

  //Camera position relative to the pixel location
  vec3 eyeVector = - position;
  vec3 eyeDirection = normalize(eyeVector);

  //Light calculations
  float lightNormDot = dot(normal, lightDirection);
  
  /////////////////////////////
  //Shadow Mapping
  /////////////////////////////
  ShadowCoord = ShadowMatrix * vec4(position, 1.0);
  float ShadowCoordW = ShadowCoord.w;
  ShadowCoord = ShadowCoord * (1.0 / ShadowCoord.w);

 //If shadow mapping is off, we want everything to be unshadowed
  float shadow = 1.0 - float(ShadowMapping);
  vec2 circle = (ShadowCoord.xy) - vec2(0.5, 0.5);

  if (bool(ShadowMapping)
      && (dot(circle, circle) < 0.25)
      && (ShadowCoord.w > 0.0))
    shadow = chebyshevUpperBound(ShadowCoord.z);

  shadow = min(shadow, smoothstep(-0.1, 1.0, lightNormDot));

  /////////////////////////////
  //Blinn Phong lighting calculation
  /////////////////////////////
  vec3 color = texture2D(colorTex, screenCoord).rgb;
  
  /////////////////////Ambient light
  float intensity = 0.2;

  vec3 ReflectedRay = reflect(-lightDirection, normal);
  intensity += 0.0001 * float(lightNormDot > 0.0) 
    * shadow * pow(max(dot(ReflectedRay, eyeDirection), 0.0), 96.0);
  
  //Scale the shadow by the shadow intensity
  shadow = 1.0 - ShadowIntensity * (1.0 - shadow);

  /////////////////////Diffuse light "shadowing"
  //The diffuse light is calculated as a "shadow", 
  float diffuseFactor = 0.5 + 0.5 * lightNormDot;
  intensity += shadow * diffuseFactor * 0.8;

  //Light attenuation
  float attenuation = min(1.0, 1.0 / (0.2 + lightDistance * (0.1 + 0.01 * lightDistance)));
  intensity *= attenuation;

  color_out = vec4(intensity * color.rgb, 1.0);
});
	}
      };
    }
  }
}

#undef STRINGIFY
