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
      }

      inline void attach(GLfloat FocalLength, GLint width, GLint height, Vector Origin)
      {
	glUseProgramObjectARB(_shaderID);
	glUniform1fARB(_FocalLengthUniform, FocalLength);
	glUniform2fARB(_WindowSizeUniform, width, height);
	glUniform3fARB(_RayOriginUniform, Origin[0], Origin[1], Origin[2]);
      }

      static inline std::string vertexShaderSource();
      static inline std::string fragmentShaderSource();
      
    protected:
      GLuint _FocalLengthUniform;
      GLuint _WindowSizeUniform;
      GLuint _RayOriginUniform;
      GLuint _depthTexUniform;
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

void main()
{
  rayDirection.x = 2.0 * gl_FragCoord.x / WindowSize.x - 1.0;
  rayDirection.y = 2.0 * gl_FragCoord.y / WindowSize.y - 1.0;
  rayDirection.y *= WindowSize.y / WindowSize.x;
  rayDirection.z = -FocalLength;
  rayDirection = (vec4(rayDirection, 0.0) * gl_ModelViewMatrix).xyz;
  rayDirection = normalize(rayDirection);

  //Cube ray intersection test
  vec3 invR = 1.0 / normalize(rayDirection);
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

  //Check what the screen depth is to make sure we don't sample into any objects
  float depth = texture2D(DepthTexture, gl_FragCoord.xy / WindowSize.xy).r;
  
  gl_FragColor = vec4(1, 0, 0.0, (tfar-tnear)/3.4642);
}
);
    }
  }
}

#include <memory>

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
  };
}
