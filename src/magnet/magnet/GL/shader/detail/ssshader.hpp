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
      namespace detail {
	/*! \brief A base class for OpenGL Shaders implementing a
	  Screen Space Filter.
	 
	  Screen space filters are filters that take a rendered image
	  and apply an image transform using only the rendered image
	  data (such as the pixel color and depth). This base class
	  simplifies the task of generating a fragment shader for each
	  pixel of the destination image.
	 */
	class SSShader : public Shader
	{
	public:
	  ~SSShader() { deinit(); }

	  /*! \brief Run the fragment shader for each pixel in the
	    output image.
	   */
	  void invoke()
	  {
	    if (!_built)
	      M_throw() << "Cannot invoke a SS filter without it being built first";

	    glDrawArrays(GL_POINTS, 0, 1);
	  }
	
	  /*!\brief An empty vertex shader. 
	    
	    All work is carried out in the Geometry shader.
	   */
	  virtual std::string initVertexShaderSource()
	  { return "#version 330\n void main() {}"; }

	virtual std::string initGeometryShaderSource()
	{
	  return
	    "#version 330\n"
	    STRINGIFY(
layout(points) in;
layout(triangle_strip) out;
layout(max_vertices = 3) out;

smooth out vec2 screenCoord;

void main()
{
  /*Here we draw a fullscreen triangle and allow the GPU to scissor to
    the screen. This prevents the difficult interpolation of the
    vertex property (screenCoord) on the diagonal of a fullscreen
    quad. This is a ridiculous optimisation I know. */

  screenCoord = vec2(0.0, 0.0);
  gl_Position = vec4(-1.0, -1.0, 0.5, 1.0);
  EmitVertex();
  
  screenCoord = vec2(2.0, 0.0);
  gl_Position = vec4(+3.0, -1.0, 0.5, 1.0);
  EmitVertex();

  screenCoord = vec2(0.0, 2.0);
  gl_Position = vec4(-1.0, +3.0, 0.5, 1.0);
  EmitVertex();
  EndPrimitive();
});
	}
	};
      }
    }
  }
}

#undef STRINGIFY
