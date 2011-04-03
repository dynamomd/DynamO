/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

    This program is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    version 3 as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#pragma once

#include "Quads.hpp"

#include <magnet/GL/detail/shader.hpp>
#include <magnet/math/vector.hpp>

namespace magnet {
  namespace GL {

    class volumeRenderer: public detail::shader<volumeRenderer>
    {
    public:      
      inline void build()
      {
	//First, call the build function in the shader
	detail::shader<volumeRenderer>::build();

	_FocalLengthUniform = glGetUniformLocationARB(_shaderID,"FocalLength");
	_WindowSizeUniform = glGetUniformLocationARB(_shaderID,"WindowSize");
	_RayOriginUniform = glGetUniformLocationARB(_shaderID,"RayOrigin");
	_depthTexUniform = glGetUniformLocationARB(_shaderID,"DepthTexture");
	_nearUniform = glGetUniformLocationARB(_shaderID,"NearDist");
	_farUniform = glGetUniformLocationARB(_shaderID,"FarDist");
	_dataTexUniform = glGetUniformLocationARB(_shaderID,"DataTexture");
      }

      inline void attach(GLfloat FocalLength, GLint width, GLint height, Vector Origin, 
			 GLint depthTex, GLint dataTex,
			 GLfloat NearDist, GLfloat FarDist)
      {
	glUseProgramObjectARB(_shaderID);
	glUniform1fARB(_FocalLengthUniform, FocalLength);
	glUniform2fARB(_WindowSizeUniform, width, height);
	glUniform3fARB(_RayOriginUniform, Origin[0], Origin[1], Origin[2]);
	glUniform1iARB(_depthTexUniform, depthTex);
	glUniform1fARB(_farUniform, FarDist);
	glUniform1fARB(_nearUniform, NearDist);
	glUniform1iARB(_dataTexUniform, dataTex);
      }

      static inline std::string vertexShaderSource();
      static inline std::string fragmentShaderSource();
      
    protected:
      GLuint _FocalLengthUniform;
      GLuint _WindowSizeUniform;
      GLuint _RayOriginUniform;
      GLuint _depthTexUniform;
      GLuint _nearUniform;
      GLuint _farUniform;
      GLuint _dataTexUniform;
    };
  }
}

//Simple macro to convert a token to a string
#define STRINGIFY(A) #A

namespace magnet {
  namespace GL {
    inline std::string 
    volumeRenderer::vertexShaderSource()
    {
      return
STRINGIFY( 
void main()
{
  gl_Position = ftransform();
}
);
    }

    inline std::string 
    volumeRenderer::fragmentShaderSource()
    {
      return
STRINGIFY( 
uniform float FocalLength;
uniform vec2 WindowSize;
uniform vec3 RayOrigin;
uniform sampler2D DepthTexture;
uniform sampler3D DataTexture;
uniform float NearDist;
uniform float FarDist;

float recalcZCoord(float zoverw)
{
  return (2.0 * NearDist * FarDist) 
    / (FarDist + NearDist - (2.0 * zoverw - 1.0) * (FarDist - NearDist));
}

void main()
{
  vec3 rayDirection;
  rayDirection.x = 2.0 * gl_FragCoord.x / WindowSize.x - 1.0;
  rayDirection.y = 2.0 * gl_FragCoord.y / WindowSize.y - 1.0;
  rayDirection.y *= WindowSize.y / WindowSize.x;
  rayDirection.z = -FocalLength;
  rayDirection = (vec4(rayDirection, 0.0) * gl_ModelViewMatrix).xyz;
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

  //If we're penetrating the volume, make sure to only cast the ray
  //from the eye position, not behind it
  if (tnear < 0) tnear = 0.0;

  vec3 tmax = max(ttop, tbot); //Distant planes
  t = min(tmax.xx, tmax.yz);//Find the first plane to be exited
  float tfar = min(t.x, t.y);//..

  //Check what the screen depth is to make sure we don't sample the volume past any objects
  float depth = recalcZCoord(texture2D(DepthTexture, gl_FragCoord.xy / WindowSize.xy).r);
  if (tfar > depth) tfar = depth;
  
  const float stepSize = 0.01;
  const float baseStepSize = 0.01;

  vec3 rayPos = RayOrigin + rayDirection * tnear;
  float length = tfar - tnear;
  vec4 color = vec4(0.0, 0.0, 0.0, 0.0);

  vec3 lightPos = (gl_ModelViewMatrixInverse * gl_LightSource[0].position).xyz;
  vec3 Specular = (gl_FrontMaterial.specular * gl_LightSource[0].specular).xyz;

  for (; (length > 0.0) && (color.a < 0.95); 
       length -= stepSize, rayPos.xyz += rayDirection * stepSize)
    {
      vec4 sample = texture3D(DataTexture, vec4((rayPos + 1.0) * 0.5, 0.0));

      vec4 src = vec4(1.0,1.0,1.0,sample.a);

      //This corrects the transparency change caused by changing step
      //size
      src.a = 1.0 - pow((1-src.a), stepSize / baseStepSize);

      //This term is to make the object more transparent
      //src.a *= 0.5;
      
      ////////////Lighting
      vec3 norm = normalize(sample.xyz * 2.0 - 1.0);
      vec3 lightDir = normalize(lightPos - rayPos);
      float lightNormDot = dot(norm, lightDir);

      //Diffuse lighting term
      float diffTerm = 0.5 * lightNormDot  + 0.5;

      //Specular lighting term
      vec3 ReflectedRay = normalize(reflect(normalize(lightDir), normalize(norm)));
      vec3 spec = 
	(lightNormDot > 0) //Test to ensure that specular is only
			   //applied to front facing voxels
	* Specular * pow(max(dot(ReflectedRay, normalize(rayDirection)), 0), 
		       gl_FrontMaterial.shininess);
      //Sum of terms 
      src.rgb = spec + diffTerm * diffTerm * src.rgb + gl_LightModel.ambient * src.rgb;
      
      
      
      ///////////Front to back blending
      src.rgb *= src.a;
      color = (1.0 - color.a) * src + color;
    }

  gl_FragColor = color;
}
);
    }
  }
}

#include <memory>
#include <magnet/GL/texture3D.hpp>

namespace coil {
  class RVolume : public RQuads
  {
  public:
    RVolume(std::string name);
    ~RVolume();
  
    virtual void initOpenGL();
    virtual void initOpenCL();
    virtual void glRender(magnet::GL::FBO& fbo);

    virtual void resize(size_t width, size_t height);

    virtual void releaseCLGLResources();
  protected:

    magnet::GL::volumeRenderer _shader;
    std::auto_ptr<magnet::GL::FBO> _fbo;

    magnet::GL::Texture3D _data;
  };
}
