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
////Quaternion mathmatics
//https://mollyrocket.com/forums/viewtopic.php?p=6154
uniform mat4 ShadowMatrix;
uniform vec3 lightPosition;

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
  normal = normalize(gl_NormalMatrix * qrot(iOrientation, vNormal.xyz));
  
  vec4 vVertex = gl_ModelViewMatrix * vec4(qrot(iOrientation, vPosition.xyz * iScale.xyz) + iOrigin.xyz, 1.0);
  gl_Position = gl_ProjectionMatrix * vVertex;
  
  ShadowCoord = ShadowMatrix * vVertex;
  lightDir =  (gl_ModelViewMatrix * vec4(lightPosition,1)).xyz - vVertex.xyz;
  eyeVector = -vVertex.xyz;
});
	}
	
	virtual std::string initFragmentShaderSource()
	{
	  return STRINGIFY(
uniform sampler2DShadow ShadowMap; //The sampler for the shadow map
uniform int ShadowMapping; //If shadow mapping is enabled or not
uniform float ShadowIntensity; //How dark the shadow is
uniform float xPixelOffset;
uniform float yPixelOffset;

varying vec4 ShadowCoord; // Texture coordinate used for shadow lookup
varying vec3 lightDir; //Direction of the light
varying vec3 normal; //The surface normal
varying vec4 color;
varying vec3 eyeVector;


void main()
{
  ///Shadow map calculations and PCF.
  float shadow = 0.0 ;
  const int steps = 3;
  const float stepsf = 3.0;
  const float stepoffset = (stepsf - 1.0) * 0.5;
  
  vec3 renormLightDir = normalize(lightDir);

  //Flip the normal for 2 sided rendering and renormalize the normals again
  vec3 renormal = normalize(normal);  
  if (!gl_FrontFacing) renormal = -renormal;
  
  //Now get the diffuse term to smooth shadow acne on back faces
  float lightNormDot = dot(renormal, renormLightDir);

  if (ShadowMapping == 1)
    {
      if (ShadowCoord.w > 0.0)
	for (int x = 0; x < steps; ++x)
	  for (int y = 0; y < steps; ++y)
	    {
	      vec4 sampleCoords = ShadowCoord;
	      sampleCoords.x 
		+= (float(x) - stepoffset) * xPixelOffset * ShadowCoord.w;
	      sampleCoords.y 
		+= (float(y) - stepoffset) * yPixelOffset * ShadowCoord.w;
	      vec2 circle = (sampleCoords.xy) / sampleCoords.w
		- vec2(0.5, 0.5);
	      
	      //Makes the light source a spotlight (and fixes a wierd
	      //line perpendicular to the light beam).
	      if (dot(circle, circle) < 0.25) 
		shadow += shadow2DProj(ShadowMap, sampleCoords).r;
	    }
	      
      shadow /= stepsf * stepsf;

      float shadeFactor = max(0.0, min(1.0, 1.0 + 2.0 * lightNormDot));
      
      //Now set a minimum on the shadow term depending on if the
      //surface is pointing away from the light source
      shadow = min(shadow, shadeFactor);
    }
  else
    shadow = 1.0;

  //Start the color calculation with the global ambient
  vec4 fragcolor = vec4(0,0,0,0);

  float lightDist = length(lightDir);

  //This is the attenuation of the light source with distance
  float attenuation = 1.0 / (1.0 + lightDist * (0.0 + 0.2 * lightDist));
 
  //Shadow intensity
  float scaledShadow = (1.0 + ShadowIntensity * (shadow - 1.0));

//  //Specular and ambient light calculation
  if (lightNormDot > 0.0)
    {
      //We use a shadow map instead of a spotlight calculation
      vec3 Eye = normalize(eyeVector);
      vec3 ReflectedRay = reflect(-renormLightDir, renormal);

      //gl_LightSource[0].specular       
      vec4 Specular = vec4(1.0,1.0,1.0,1.0) * pow(max(dot(ReflectedRay, Eye), 0.0), 25);
      //We multiply by shadow to stop the specular highlight appearing in the shadow
      fragcolor += shadow * attenuation * Specular;
    }

  fragcolor += attenuation * vec4(0.2);

  //The diffusive term
  float rescaledDot = 0.5 * lightNormDot + 0.5;
  fragcolor += attenuation * scaledShadow * rescaledDot * rescaledDot * color;

  fragcolor.a = color.a;

  gl_FragColor = fragcolor;
});
	}
      };
    }
  }
}

#undef STRINGIFY
