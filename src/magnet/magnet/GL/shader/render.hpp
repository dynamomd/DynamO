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
#include <iostream>
#define STRINGIFY(A) #A

namespace magnet {
  namespace GL {
    namespace shader {
      /*! \brief An instancing shadow mapping/diffuse/ambient and specular shader
       * with per-pixel lighting for rendering scenes.
       */
      class RenderShader: public detail::Shader
      {
      public:
	virtual std::string initVertexShaderSource()
	{
	  return STRINGIFY(
uniform mat4 ShadowMatrix;
uniform vec3 lightPosition;

uniform mat4 ProjectionMatrix;
uniform mat4 ViewMatrix;
uniform mat3 NormalMatrix;

attribute vec4 vPosition;
attribute vec4 vColor;
attribute vec4 vNormal;
attribute vec4 iOrigin;
attribute vec4 iOrientation;
attribute vec4 iScale;

varying vec4 ShadowCoord;
varying vec3 lightDir;
varying vec3 normal;
varying vec4 color;
varying vec3 eyeVector;

vec3 qrot(vec4 q, vec3 v)
{ return v + 2.0 * cross(cross(v,q.xyz) + q.w * v, q.xyz); } 

void main()
{
  color = vColor;
  normal = normalize(NormalMatrix * qrot(iOrientation, vNormal.xyz));
  
  vec3 scale = iScale.xyz + vec3(equal(iScale.xyz, vec3(0.0))) * iScale.x;

  vec4 vVertex = ViewMatrix * vec4(qrot(iOrientation, vPosition.xyz * scale) + iOrigin.xyz, 1.0);
  gl_Position = ProjectionMatrix * vVertex;
  
  ShadowCoord = ShadowMatrix * vVertex;
  lightDir =  (ViewMatrix * vec4(lightPosition,1)).xyz - vVertex.xyz;
  eyeVector = -vVertex.xyz;
});
	}
	
	virtual std::string initFragmentShaderSource()
	{
	  return STRINGIFY(
uniform sampler2D ShadowMap; //The sampler for the shadow map
uniform int ShadowMapping; //If shadow mapping is enabled or not
uniform float ShadowIntensity; //How dark the shadow is

varying vec4 ShadowCoord; // Texture coordinate used for shadow lookup
varying vec3 lightDir; //Direction of the light
varying vec3 normal; //The surface normal
varying vec4 color;
varying vec3 eyeVector;

vec4 ShadowCoordPostW;

float chebyshevUpperBound(float distance)
{
  vec2 moments = texture2D(ShadowMap,ShadowCoordPostW.xy).rg;
	
  // Surface is fully lit. as the current fragment is before the light occluder
  if (distance <= moments.x)
    return 1.0 ;

  // The fragment is either in shadow or penumbra. We now use
  // chebyshev's upperBound to check How likely this pixel is to be
  // lit (p_max)
  float variance = moments.y - (moments.x * moments.x);
  variance = max(variance, 0.0000001);

  float d = distance - moments.x;
  float p_max = variance / (variance + d * d);

  return p_max;
}

void main()
{
  ///Shadow map calculations and PCF.
  const int steps = 3;
  const float stepsf = float(steps);
  const float stepoffset = (stepsf - 1.0) * 0.5;
  
  //Flip the normal for 2 sided rendering and renormalize the normals again
  vec3 renormal = normalize(normal);  
  if (!gl_FrontFacing) renormal = -renormal;
  
  //Now get the diffuse term to smooth shadow acne on back faces
  vec3 renormLightDir = normalize(lightDir);
  float lightNormDot = dot(renormal, renormLightDir);

  ShadowCoordPostW = ShadowCoord / ShadowCoord.w;

  //If shadow mapping is off, we want everything to be unshadowed
  float shadow = 1 - bool(ShadowMapping);
  vec2 circle = (ShadowCoordPostW.xy) - vec2(0.5, 0.5);

  if (bool(ShadowMapping)
      && (dot(circle, circle) < 0.25)
      && (ShadowCoord.w > 0.0))
    shadow = chebyshevUpperBound(ShadowCoordPostW.z);
  
  //This term accounts for self shadowing, to help remove shadow acne
  //on the sides of objects
  shadow = min(shadow, smoothstep(-0.1, 1.0, lightNormDot));

  /////////////////////Ambient light
  float intensity = 0.2;

  /////////////////////Specular light
  //This is calculated before scaling the shadow, as specular lighting
  //should never appear in the shadow.
  vec3 Eye = normalize(eyeVector);
  vec3 ReflectedRay = reflect(-renormLightDir, renormal);
  intensity += float(lightNormDot > 0.0) * shadow * pow(max(dot(ReflectedRay, Eye), 0.0), 96.0);
  
  //Scale the shadow by the shadow intensity
  shadow = 1.0 - ShadowIntensity * (1.0 - shadow);

  /////////////////////Diffuse light "shadowing"
  //The diffuse light is calculated as a "shadow", 
  float diffuseFactor = 0.5 + 0.5 * lightNormDot;
  intensity += shadow * diffuseFactor * 0.8;

  //Light attenuation
  float lightDist = length(lightDir);
  float attenuation = min(1.0, 1.0 / (0.2 + lightDist * (0.1 + 0.01 * lightDist)));
  intensity *= attenuation;

  gl_FragColor = vec4(intensity * color.rgb, color.a);
});
	}
      };
    }
  }
}

#undef STRINGIFY
