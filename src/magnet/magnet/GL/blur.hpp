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

    class blurFilter : public detail::filter<blurFilter, 5>
    {
    public:
      inline static const GLfloat* weights()
      {
	static const GLfloat weights[5][5] = 
	  {
	    {1.0/331.0, 4/331.0, 7/331.0, 4/331.0, 1/331.0},
	    {4/331.0, 20/331.0, 33/331.0, 20/331.0, 4/331.0},
	    {7/331.0, 33/331.0, 55/331.0, 33/331.0, 7/331.0},
	    {4/331.0, 20/331.0, 33/331.0, 20/331.0, 4/331.0},
	    {1/331.0, 4/331.0, 7/331.0, 4/331.0, 1/331.0}
	  };
	
	return (const GLfloat*)weights;
      }
	
    };

    class boxFilter : public detail::filter<boxFilter, 5>
    {
    public:
      inline static const GLfloat* weights()
      {
	static const GLfloat weights[5][5] = 
	  {
	    {1.0/331.0, 4/331.0, 7/331.0, 4/331.0, 1/331.0},
	    {4/331.0, 20/331.0, 33/331.0, 20/331.0, 4/331.0},
	    {7/331.0, 33/331.0, 55/331.0, 33/331.0, 7/331.0},
	    {4/331.0, 20/331.0, 33/331.0, 20/331.0, 4/331.0},
	    {1/331.0, 4/331.0, 7/331.0, 4/331.0, 1/331.0}
	  };
	
	return (const GLfloat*)weights;
      }
	
    };
  }
}

