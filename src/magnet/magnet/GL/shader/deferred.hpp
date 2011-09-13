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
	  return "#version 110\n" STRINGIFY(
uniform sampler2D colorTex;
uniform sampler2D normalTex;
uniform sampler2D positionTex;
uniform vec3 lightPosition;
uniform vec3 camPosition;
varying vec2 screenCoord;

void main()
{
  vec3 color = texture2D(colorTex, screenCoord).rgb;
  vec3 normal = texture2D(normalTex, screenCoord).rgb;
  normal = normalize(normal);

  //Model space position of the vertex
  vec3 position = texture2D(positionTex, screenCoord).xyz;
  
  //light position relative to the pixel location
  vec3 lightVector = lightPosition - position;
  float lightDistance = length(lightVector);
  vec3 lightDirection = lightVector * (1.0 / lightDistance);

  //Camera position relative to the pixel location
  vec3 eyeVector = camPosition - position;
  vec3 eyeDirection = normalize(eyeVector);

  //Light calculations
  float lightNormDot = dot(normal, lightDirection);

  float shadow = 1.0;

  /////////////////////Ambient light
  float intensity = 0.2;

  vec3 ReflectedRay = reflect(-lightDirection, normal);
  intensity += 0.0001 * float(lightNormDot > 0.0) * shadow * pow(max(dot(ReflectedRay, eyeDirection), 0.0), 96.0);
  
  //Scale the shadow by the shadow intensity
  //shadow = 1.0 - ShadowIntensity * (1.0 - shadow);

  /////////////////////Diffuse light "shadowing"
  //The diffuse light is calculated as a "shadow", 
  float diffuseFactor = 0.5 + 0.5 * lightNormDot;
  intensity += shadow * diffuseFactor * 0.8;

  //Light attenuation
  float attenuation = min(1.0, 1.0 / (0.2 + lightDistance * (0.1 + 0.01 * lightDistance)));
  intensity *= attenuation;

  gl_FragData[0] = vec4(intensity * color.rgb, 1.0);
});
	}
      };
    }
  }
}

#undef STRINGIFY
