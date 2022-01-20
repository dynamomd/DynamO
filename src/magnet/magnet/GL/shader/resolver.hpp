/*    dynamo:- Event driven molecular dynamics simulator 
 *    http://www.dynamomd.org
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
      class ResolverShader: public detail::SSShader
      {
      public:
	virtual std::string initFragmentShaderSource()
	{
	  return STRINGIFY(
layout (location = 0) out vec4 outTex;
uniform sampler2DMS inTex;
uniform int sample;

void main()
{
  outTex = texelFetch(inTex, ivec2(gl_FragCoord.xy), sample);
});
	}
      };
    }
  }
}

#undef STRINGIFY
