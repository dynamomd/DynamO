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
#include <magnet/GL/shader/detail/filter.hpp>

namespace magnet {
  namespace GL {
    namespace shader {
      /*! \brief Implements a 3x3 HiPass Shader.
	
	These screen space filters sharpen images by only allowing
	high frequency data through.
       */
      class HiPass3x3 : public detail::SSKernelShader
      {
      public:
	void build() { SSKernelShader::build(3); }

	virtual const GLfloat* weights()
	{
	  static const GLfloat weights[3][3] = 
	    {
	      {-1, -1, -1},
	      {-1,  9, -1},
	      {-1, -1, -1}
	    };
	
	  return (const GLfloat*)weights;
	}
	
      };

      /*! \brief Implements a 5x5 HiPass filter. */
      class HiPass5x5 : public detail::SSKernelShader
      {
      public:
	void build() { SSKernelShader::build(5); }
      
	virtual const GLfloat* weights()
	{
	  static const GLfloat weights[5][5] = 
	    {
	      {-1/1.0, -1/1.0, -1/1.0, -1/1.0, -1/1.0},
	      {-1/1.0, -1/1.0, -1/1.0, -1/1.0, -1/1.0},
	      {-1/1.0, -1/1.0, 25/1.0, -1/1.0, -1/1.0},
	      {-1/1.0, -1/1.0, -1/1.0, -1/1.0, -1/1.0},
	      {-1/1.0, -1/1.0, -1/1.0, -1/1.0, -1/1.0}
	    };
	
	  return (const GLfloat*)weights;
	}	
      };
    }
  }
}

