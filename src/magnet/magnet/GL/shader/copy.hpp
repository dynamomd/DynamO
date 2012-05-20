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
      class CopyShader : public detail::SSKernelShader
      {
      public:
	void build() { SSKernelShader::build(1); }

	virtual const GLfloat* weights()
	{
	  static const GLfloat weights[1][1] = { {1} };
	  return (const GLfloat*)weights;
	}
      };
    }
  }
}
