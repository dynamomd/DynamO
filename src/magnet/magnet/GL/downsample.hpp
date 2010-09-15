/*    DYNAMO:- Event driven molecular dynamics simulator 
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

#include <magnet/GL/detail/filter.hpp>

namespace magnet {
  namespace GL {    
    class downsampleFilter: public detail::filter<downsampleFilter>
    {
    public:      
      void build(GLuint FBO, GLsizei width, GLsizei height)
      {
	if (!glewIsSupported("GL_ARB_half_float_pixel"))
	  M_throw() << "Floating point textures not supported!";

	detail::filter<downsampleFilter>::build(FBO, width, height, GL_RGBA16F_ARB, GL_FLOAT);
      }

      void build(GLsizei width, GLsizei height)
      {
	if (!glewIsSupported("GL_ARB_half_float_pixel"))
	  M_throw() << "Floating point textures not supported!";

	detail::filter<downsampleFilter>::build(width, height, GL_RGBA16F_ARB, GL_FLOAT);
      }
      
      static inline std::string vertexShaderSource();
      static inline std::string fragmentShaderSource();
    

      void invoke()
      {
	
      }
    protected:
      GLuint _OutputTexture;
      
    };
  }
}

#include <magnet/GL/detail/shaders/downsample.glh>
