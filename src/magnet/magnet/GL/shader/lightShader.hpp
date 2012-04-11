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
      class PointLightShader: public detail::SSShader
      {
      public:
	virtual std::string initFragmentShaderSource()
	{
	  return STRINGIFY(
layout (location = 0) out vec4 color_out;

//Standard G-buffer data
uniform sampler2DMS colorTex;
uniform sampler2DMS normalTex;
uniform sampler2DMS positionTex;
uniform vec3 lightPosition;
uniform vec3 lightColor;
uniform float ambientLight;
uniform float lightAttenuation;
uniform float lightSpecularExponent;
uniform float lightSpecularFactor;
uniform int samples;

vec3 calcLighting(vec3 position, vec3 normal, vec3 diffuseColor)
{  
  vec3 lightVector = lightPosition - position;
  float lightDistance = length(lightVector);
  vec3 lightDirection = lightVector * (1.0 / lightDistance);
 
  //if the normal has a zero length, illuminate it as though it was
  //fully lit
  float normal_length = length(normal);
  normal = (normal_length == 0) ?  lightDirection : normal / normal_length;
 
  float lightNormDot = dot(normal, lightDirection);

  /////////////////////////////
  //Blinn Phong lighting calculation
  /////////////////////////////

  vec3 ReflectedRay = reflect(-lightDirection, normal);

  vec3 eyeDirection = normalize(-position);
  //Specular
  float specular = lightSpecularFactor * float(lightNormDot > 0.0)
    * pow(max(dot(ReflectedRay, eyeDirection), 0.0), lightSpecularExponent);
  
  float diffuse = smoothstep(-0.5, 1.0, lightNormDot);

  //Light attenuation
  float decay_factor = 1.0 / (1.0 + lightAttenuation * lightDistance * lightDistance);

  return decay_factor * lightColor * (specular + diffuse * diffuseColor);
}

void main()
{
  //Now calculate the color from the samples
  vec4 color_sum = vec4(0.0);
  
  for (int sample_id = 0; sample_id < samples; sample_id++)
    {
      vec4 color = texelFetch(colorTex, ivec2(gl_FragCoord.xy), sample_id).rgba;

      //If alpha is zero, this is an empty pixel, and should not
      //contribute to the tone mapping
      if (color.a != 0)
	{
	  //Eye space normal and position
	  vec3 normal = texelFetch(normalTex, ivec2(gl_FragCoord.xy), sample_id).rgb;
	  vec3 position = texelFetch(positionTex, ivec2(gl_FragCoord.xy), sample_id).xyz;
	  color_sum.rgb += ambientLight * color.rgb + calcLighting(position, normal, color.rgb);
	  color_sum.a += 1.0;
	}
    }
 
  //We write out the HDR color here, along with the occupancy
  //(fraction of drawn pixels) in the alpha channel.
  color_out = color_sum / float(samples);
});
	}
      };
    }
  }
}

#undef STRINGIFY
