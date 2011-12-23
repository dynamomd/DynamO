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
      class HDRCombinerShader: public detail::SSShader
      {
      public:
	virtual std::string initFragmentShaderSource()
	{
	  return "#version 330\n"
	    STRINGIFY(
//Normalized position on the screen
smooth in vec2 screenCoord;
layout (location = 0) out vec4 color_out;

//Standard G-buffer data
uniform sampler2D colorTex;
uniform sampler2DMS depthTex;
uniform float invGamma;
uniform float exposure;

void main()
{
  ivec2 pixelcoord = ivec2(textureSize(depthTex) * screenCoord);

  //Copy the first sample depth across
  gl_FragDepth = texelFetch(depthTex, pixelcoord, 0).r;

  vec4 color = texelFetch(colorTex, pixelcoord, 0).rgba;

  color_out = vec4(pow(color.rgb * exposure, vec3(invGamma)), 1.0);
});
	}
      };
    }
  }
}

#undef STRINGIFY
