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
      }

      inline void attach(GLfloat FocalLength, GLint width, GLint height, Vector Origin)
      {
	glUseProgramObjectARB(_shaderID);
	glUniform1fARB(_FocalLengthUniform, FocalLength);
	glUniform2iARB(_WindowSizeUniform, width, height);
	glUniform3fARB(_RayOriginUniform, Origin[0], Origin[1], Origin[2]);
      }

      static inline std::string vertexShaderSource();
      static inline std::string fragmentShaderSource();
      
    protected:
      GLuint _FocalLengthUniform;
      GLuint _WindowSizeUniform;
      GLuint _RayOriginUniform;
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

void main()
{
  gl_FragColor = vec4(1,0,0,1);
}
);
    }
  }
}

namespace coil {
  class RVolume : public RQuads
  {
  public:
    RVolume(std::string name);
  
    virtual void initOpenGL();
    virtual void initOpenCL();
    virtual void glRender();

  protected:

    magnet::GL::volumeRenderer _shader;
  };
}
