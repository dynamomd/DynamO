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
#define STRINGIFY(A) #A

namespace magnet {
  namespace GL {
    namespace shader {
      class DOF : public detail::Shader
      {
      public:
	void invoke()
	{
	  //Setup the shader arguments
	  glUseProgram(_shaderID);
	  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	  drawScreenQuad();
	  //Restore the fixed pipeline
	  glUseProgramObjectARB(0);
	}
	
	virtual std::string fragmentShaderSource()
	{
	  return STRINGIFY(
uniform sampler2D u_Texture0; //Blurred image
uniform sampler2D u_Texture1; //Original
uniform sampler2D u_Texture2; //Depth buffer
uniform float focalDistance;
uniform float focalRange;
uniform float nearDist;
uniform float farDist;

float LinearizeDepth(float zoverw)
{
  return(2.0 * nearDist) / (farDist + nearDist - zoverw * (farDist - nearDist));
}

void main(void)
{
  float fcldist = focalDistance;
  if (fcldist == 0) //Automatic mode
    fcldist = LinearizeDepth(texture2D(u_Texture2, vec2(0.5,0.5)).r);
  
  vec4 original = texture2D(u_Texture1, gl_TexCoord[0].st);
  vec4 blurred = texture2D(u_Texture0, gl_TexCoord[0].st);
  
  float depth = LinearizeDepth(texture2D(u_Texture2, gl_TexCoord[0].st).r);
  float blur = clamp(abs(depth - fcldist) / focalRange, 0.0, 1.0);
  
  //gl_FragColor =  vec4(blur,0,0,0);
  gl_FragColor = original + blur * (blurred - original);
});
	}
      };
    }
  }
}
#undef STRINGIFY
