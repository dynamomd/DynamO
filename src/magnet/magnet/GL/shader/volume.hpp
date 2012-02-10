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

//Simple macro to convert a token to a string
#define STRINGIFY(A) #A

namespace magnet {
  namespace GL {
    namespace shader {
      /*! \brief A shader for ray-tracing cubic volumes.
       *
       * This shader will render a volume data set in one pass. For
       * more information on the method please see
       * https://www.marcusbannerman.co.uk/index.php/home/42-articles/97-vol-render-optimizations.html
       */
      class VolumeShader: public detail::Shader
      {
      public:      
	virtual std::string initVertexShaderSource()
	{ return "#version 330\n"
	    STRINGIFY(
uniform mat4 ProjectionMatrix;
uniform mat4 ViewMatrix;

layout (location = 0) in vec4 vPosition;
layout (location = 3) in vec4 iOrigin;
layout (location = 4) in vec4 iOrientation;
layout (location = 5) in vec4 iScale;

vec3 qrot(vec4 q, vec3 v)
{ return v + 2.0 * cross(cross(v,q.xyz) + q.w * v, q.xyz); } 

void main()
{ 
  vec3 scale = iScale.xyz + vec3(equal(iScale.xyz, vec3(0.0))) * iScale.x;
  vec4 vVertex
    = ViewMatrix
    * vec4(qrot(iOrientation, vPosition.xyz * scale) + iOrigin.xyz, 1.0);

  gl_Position = ProjectionMatrix * vVertex; }
); 
	}

	virtual std::string initFragmentShaderSource()
	{
	  return "#version 330\n"
	    STRINGIFY(
uniform float FocalLength;
uniform vec2 WindowSize;
uniform vec3 RayOrigin;

//The light's position in model space (untransformed)
uniform sampler1D TransferTexture;
uniform sampler1D IntTransferTexture;
uniform sampler2D DepthTexture;
uniform sampler3D DataTexture;
uniform float StepSize;
uniform float DitherRay;

uniform mat4 ProjectionMatrix;
uniform mat4 ViewMatrix;

layout (location = 0) out vec4 color_out;

//This function converts a value in a depth buffer back into a object-space distance
float recalcZCoord(float zoverw)
{
  float A = ProjectionMatrix[2].z;
  float B = ProjectionMatrix[3].z;
  float zNearDist =  -B / (1.0 - A);
  float zFarDist = B / (1.0 + A);
  
  return (2.0 * zNearDist * zFarDist) 
    / (zFarDist + zNearDist - (2.0 * zoverw - 1.0) * (zFarDist - zNearDist));
}

uniform vec3 lightPosition;
uniform vec3 lightColor;
uniform float ambientLight;
uniform float lightIntensity;
uniform float lightAttenuation;
uniform float lightSpecularExponent;
uniform float lightSpecularFactor;

vec3 calcLighting(vec3 position, vec3 normal, vec3 diffuseColor)
{  
  vec3 lightVector = lightPosition - position;
  float lightDistance = length(lightVector);
  vec3 lightDirection = lightVector * (1.0 / lightDistance);
 
  //if the normal has a zero length, illuminate it as though it was
  //not lit
  float normal_length = length(normal);
  normal = (normal_length == 0) ?  -lightDirection : normal / normal_length;
 
  //Camera position relative to the pixel location
  vec3 eyeVector = -position;
  vec3 eyeDirection = normalize(eyeVector);

  float lightNormDot = dot(normal, lightDirection);
  //Light attenuation
  float intensity = lightIntensity / (1.0 + lightAttenuation * lightDistance * lightDistance);

  /////////////////////////////
  //Blinn Phong lighting calculation
  /////////////////////////////

  vec3 ReflectedRay = reflect(-lightDirection, normal);
  //Specular
  float specular = lightSpecularFactor * float(lightNormDot > 0.0)
    * pow(max(dot(ReflectedRay, eyeDirection), 0.0), lightSpecularExponent);
  
  float diffuse = smoothstep(-1.1, 1.0, lightNormDot);

  return intensity 
    * (specular * lightColor
       + diffuse * diffuseColor * lightColor);
}


void main()
{
  //Calculate the ray direction using viewport information
  vec3 rayDirection;
  rayDirection.x = 2.0 * gl_FragCoord.x / WindowSize.x - 1.0;
  rayDirection.y = 2.0 * gl_FragCoord.y / WindowSize.y - 1.0;
  rayDirection.y *= WindowSize.y / WindowSize.x;
  rayDirection.z = -FocalLength;
  rayDirection = (vec4(rayDirection, 0.0) * ViewMatrix).xyz;
  rayDirection = normalize(rayDirection);
  
  //Cube ray intersection test
  vec3 invR = 1.0 / rayDirection;
  vec3 boxMin = vec3(-1.0,-1.0,-1.0);
  vec3 boxMax = vec3( 1.0, 1.0, 1.0);
  vec3 tbot = invR * (boxMin - RayOrigin);
  vec3 ttop = invR * (boxMax - RayOrigin);
  
  //Now sort all elements of tbot and ttop to find the two min and max elements
  vec3 tmin = min(ttop, tbot); //Closest planes
  vec2 t = max(tmin.xx, tmin.yz); //Out of the closest planes, find the last to be entered (collision point)
  float tnear = max(t.x, t.y);//...
  
  //If the viewpoint is penetrating the volume, make sure to only cast the ray
  //from the eye position, not behind it
  if (tnear < 0.0) tnear = 0.0;
  
  //Now work out when the ray will leave the volume
  vec3 tmax = max(ttop, tbot); //Distant planes
  t = min(tmax.xx, tmax.yz);//Find the first plane to be exited
  float tfar = min(t.x, t.y);//...
  
  //Check what the screen depth is to make sure we don't sample the
  //volume past any standard GL objects
  float bufferDepth = texture(DepthTexture, gl_FragCoord.xy / WindowSize.xy).r;
  float depth = recalcZCoord(bufferDepth);
  if (tfar > depth) tfar = depth;
  
  //We need to calculate the ray's starting position. We add a random
  //fraction of the stepsize to the original starting point to dither
  //the output
  float random = DitherRay * fract(sin(gl_FragCoord.x * 12.9898 + gl_FragCoord.y * 78.233) * 43758.5453); 
  vec3 rayPos = RayOrigin + rayDirection * (tnear + StepSize * random);
  
  //The color accumulation variable
  vec4 color = vec4(0.0, 0.0, 0.0, 0.0);
  
  //Start the sampling by initialising the ray variables
  float lastsamplea = 0.0;
  vec4 lastTransfer = texture(IntTransferTexture, lastsamplea);
  vec3 lastnorm = vec3(0.0, 0.0, 0.0); 
  
  for (float length = tfar - tnear; length > 0.0; 
       length -= StepSize, rayPos.xyz += rayDirection * StepSize)
    {
      //Grab the volume sample
      vec4 sample = texture(DataTexture, (rayPos + 1.0) * 0.5);
  
//      //Sort out the normal data
      vec3 norm = sample.xyz * 2.0 - 1.0;
      norm = normalize(rayPos);
      //Test if we've got a bad normal and need to reuse the old one
      float sqrnormlength = dot(norm,norm);
      norm /= sqrt(sqrnormlength);
      if (sqrnormlength < 0.01)
	norm = lastnorm; 

      //Store the current normal
      lastnorm = norm;      	
      norm = (ViewMatrix * vec4(norm, 0.0)).xyz;

      vec4 src = vec4(0.0, 0.0, 0.0, 0.0);
      float delta = sample.a - lastsamplea;
      vec4 transfer = texture(IntTransferTexture, sample.a);
      float deltaT = transfer.a - lastTransfer.a;
      vec3 deltaK = transfer.rgb - lastTransfer.rgb;

      if (delta == 0.0)
	{ //Special case where the integration breaks down, just use the constant val.
	  src = texture(TransferTexture, sample.a);
	  src.a = (1.0 - exp( - StepSize * src.a));
	}
      else
	{
	  /*Pre-Integrated color calc*/
	  float opacity = 1.0 - exp( - deltaT * StepSize / delta);	  
	  vec3 color = abs(deltaK) / (abs(deltaT) + 1.0e-10);
	  src = vec4(color, opacity);
	}

      lastTransfer = transfer;
      lastsamplea = sample.a;

      ////////////Lighting calculations
      //We perform all the calculations in the eye space
      src.rgb = calcLighting((ViewMatrix * vec4(rayPos,1.0)).xyz, norm, src.rgb);
      ///////////Front to back blending
      src.rgb *= src.a;
      color = (1.0 - color.a) * src + color;
  
      //We only accumulate up to 0.95 alpha (the blending never
      //reaches 1).
      if (color.a >= 0.95)
  	{
  	  //We have to renormalize the color by the alpha value (see
  	  //below)
  	  color.rgb /= color.a;
  	  //Set the alpha to one to make sure the pixel is not transparent
  	  color.a = 1.0;
  	  break;
  	}
    }
  /*We must renormalize the color by the alpha value. For example, if
  our ray only hits just one white voxel with a alpha of 0.5, we will have
  
  src.rgb = vec4(1,1,1,0.5)
  
  src.rgb *= src.a; 
  //which gives, src.rgb = 0.5 * src.rgb = vec4(0.5,0.5,0.5,0.5)
  
  color = (1.0 - color.a) * src + color;
  //which gives, color = (1.0 - 0) * vec4(0.5,0.5,0.5,0.5) + vec4(0,0,0,0) = vec4(0.5,0.5,0.5,0.5)
  
  So the final color of the ray is half way between white and black, but the voxel it hit was white!
  The solution is to divide by the alpha, as this is the "amount of color" added to color.
  */
  color.rgb /= float(color.a == 0.0) + color.a;
  color_out = color;
});
	}
      };
    }
  }
}

#undef STRINGIFY
