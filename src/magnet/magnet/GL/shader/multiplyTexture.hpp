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
      /*! \brief Simple \ref Shader to multiply two textures.*/
      class MultiplyTexture : public detail::SSShader
      {
      public:
	virtual std::string initFragmentShaderSource()
	{
	  return "#version 330\n"
	    STRINGIFY(
uniform sampler2D u_Texture0; //input
uniform sampler2D u_Texture1; //Depth buffer

smooth in vec2 screenCoord;
layout (location = 0) out vec4 color_out;

void main(void)
{
  color_out 
    = texture(u_Texture0, screenCoord) 
    * texture(u_Texture1, screenCoord);
});
	}
      };
    }
  }
}

#undef STRINGIFY
