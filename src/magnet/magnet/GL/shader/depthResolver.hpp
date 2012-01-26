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
      class DepthResolverShader: public detail::SSShader
      {
      public:
	virtual std::string initFragmentShaderSource()
	{
	  return "#version 330\n"
	    STRINGIFY(
layout (location = 0) out vec4 color_out;
uniform sampler2DMS posTex;
uniform int samples;
uniform mat4 ViewMatrix;
uniform mat4 ProjectionMatrix;

void main()
{
  //Now calculate the color from the samples
  float out_depth = 1.0;
  for (int sample_id = 0; sample_id < samples; sample_id++)
    {
      vec4 pos = vec4(texelFetch(posTex, ivec2(gl_FragCoord.xy), sample_id).xyz,
		      1.0);
      pos = ViewMatrix * pos;
      pos = ProjectionMatrix * pos;
      pos = pos / pos.w;
      
      out_depth = min((pos.z + 1.0) / 2.0, out_depth);
    }
  
  gl_FragDepth = out_depth;

  color_out = vec4(1.0);
});
	}
      };
    }
  }
}

#undef STRINGIFY
