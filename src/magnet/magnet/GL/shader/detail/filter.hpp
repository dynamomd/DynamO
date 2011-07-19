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
#include <sstream>
#define STRINGIFY(A) #A

namespace magnet {
  namespace GL {
    namespace shader {
      namespace detail {
	/*! \brief A base class for OpenGL Shaders implementing a
	 * Screen Space filter using a square kernel.
	 *
	 * This class can be used to implement simple kernel-based
	 * effects like a Gaussian blur, box filter or laplacian.
	 *
	 * Most simple screen space filters typically take the form of
	 * a "kernel". This kernel takes a square of pixels
	 * surrounding the input pixel, and \ref weights() these
	 * surrounding pixels together to calculate the output pixel.
	 *
	 * Derived classes only need to provide the size
	 * of the kernel to the constructor, and define the \ref
	 * weights() function which specifies a square array
	 * specifying the scaling factors applied to the surrounding
	 * pixels before they are summed together to form the output
	 * pixel.
	 */
	class SSKernelShader : public SSShader
	{
	protected:
	  /*! \brief Builds the screen space shader and allocates the
	   * associated OpenGL objects.
	   *
	   * This calls the underlying \ref Shader::build() function,
	   * and then loads the weights uniform with the kernel of the
	   * filter.
	   *
	   * \param stencilwidth The width/height of the filter kernel.
	   */
	  void build(int stencilwidth)
	  {
	    _stencilwidth = stencilwidth;

	    Shader::build();
	    //Get the shader args
	    glUseProgram(_shaderID);
	    //Set the weights now
	    GLint weightsUniform = glGetUniformLocationARB(_shaderID, "weights");
	    glUniform1fvARB(weightsUniform, stencilwidth * stencilwidth, weights());
	  }

	public:
	  /*! \brief Generates the fragment shader source according to
	   * the kernel's size.
	   */
	  virtual std::string initFragmentShaderSource()
	  {
	    //Simple writethrough fragment shader
	    std::ostringstream data;
	    data << _stencilwidth; 
	  
	    return std::string("#define stencilwidth ") + data.str() + "\n"
STRINGIFY(
uniform vec2 u_Scale;
uniform float weights[stencilwidth * stencilwidth];
uniform sampler2D u_Texture0;

varying vec2 screenCoord;

void main()
{
  vec4 color = vec4(0.0, 0.0, 0.0, 0.0);
  for(int x = 0; x < stencilwidth; ++x)
    for(int y = 0; y < stencilwidth; ++y)
      color += weights[y * stencilwidth + x] * texture2D(u_Texture0, screenCoord
							 + vec2((x - stencilwidth / 2) * u_Scale.x, 
								(y - stencilwidth / 2) * u_Scale.y));
  
  gl_FragColor = vec4(color.rgb, 1.0);
});
	  }

	protected:
	  /*! \brief Returns the weights that make up the kernel of
	   * the filter.
	   *
	   * This needs to be a _stencilwidth x _stencilwidth array of
	   * floats and the sum of the elements of the array should
	   * equal 1 (for a normalized filter).
	   */
	  virtual const GLfloat* weights() = 0;

	  int _stencilwidth; 
	}; 
      }
    } 
  } 
}

#undef STRINGIFY
