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

    class laplacianFilter : public detail::filter<laplacianFilter, 5>
    {
    public:
      inline static const GLfloat* weights()
      {
	static const GLfloat weights[5][5] = 
	  {
	    { 0,  0, -1,  0,  0},
	    { 0, -1, -2, -1,  0},
	    {-1, -2, 16, -2, -1},
	    { 0, -1, -2, -1,  0},
	    { 0,  0, -1,  0,  0}
	  };
	
	return (const GLfloat*)weights;
      }
	
    };
  }
}
