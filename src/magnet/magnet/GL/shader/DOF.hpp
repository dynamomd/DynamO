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
      /*! \brief A Depth Of Field Shader.
       *
       * This shader will give a depth of field effect, by blending
       * two textures together according to the pixel depth.
       *
       * The u_Texture0 uniform should contain a very "blurred" image.
       *
       * The u_Texture1 uniform should contain a "sharp" image.
       *
       * The u_Texture2 uniform should contain the depth information of the scene.
       *
       * For each pixel, the depth is looked up in u_Texture2. If this
       * pixel is in focus (set by the focalDistance uniform) then it
       * will be sampled from u_Texture1. If the pixel is out of focus
       * it is sampled from u_Texture0. The two textures are smoothly
       * blended together over a range set by the focalRange uniform.
       */
      class DOFShader : public detail::SSShader
      {
      public:
	/*! \brief The actual DOF filter. */
	virtual std::string initFragmentShaderSource()
	{
	  return "#version 330\n"
	    STRINGIFY(
uniform sampler2D u_Texture0; //Blurred image
uniform sampler2D u_Texture1; //Original
uniform sampler2D u_Texture2; //Depth buffer
uniform float focalDistance;
uniform float focalRange;
uniform float nearDist;
uniform float farDist;

smooth in vec2 screenCoord;
layout (location = 0) out vec4 color_out;

float LinearizeDepth(float zoverw)
{
  return(2.0 * nearDist) / (farDist + nearDist - zoverw * (farDist - nearDist));
}

void main(void)
{
  float fcldist = focalDistance;
  if (fcldist == 0) //Automatic mode
    fcldist = LinearizeDepth(texture(u_Texture2, vec2(0.5,0.5)).r);
  
  vec4 original = texture(u_Texture1, screenCoord);
  vec4 blurred = texture(u_Texture0, screenCoord);
  
  float depth = LinearizeDepth(texture(u_Texture2, screenCoord).r);
  float blur = clamp(abs(depth - fcldist) / focalRange, 0.0, 1.0);
  
  color_out = original + blur * (blurred - original);
});
	}
      };
    }
  }
}
#undef STRINGIFY
