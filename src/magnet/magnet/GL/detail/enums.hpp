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

namespace magnet {
  namespace GL {
    /*! \brief Enumerations describing the type of elements that can
     * be rendered.
     */
    namespace element_type {
      enum Enum 
	{
	  POINTS = GL_POINTS,
	  LINE_STRIP = GL_LINE_STRIP,
	  LINE_LOOP = GL_LINE_LOOP,
	  LINES = GL_LINES,
	  TRIANGLE_STRIP = GL_TRIANGLE_STRIP,
	  TRIANGLE_FAN = GL_TRIANGLE_FAN,
	  TRIANGLES = GL_TRIANGLES,
	  QUAD_STRIP = GL_QUAD_STRIP,
	  QUADS = GL_QUADS,
	  POLYGON = GL_POLYGON
	};
    }

    //! \brief The available GL targets to which this Buffer may be
    //! bound.
    namespace buffer_targets {
      enum Enum
	{ ARRAY = GL_ARRAY_BUFFER,
	  ELEMENT_ARRAY = GL_ELEMENT_ARRAY_BUFFER, 
	  PIXEL_PACK_BUFFER = GL_PIXEL_PACK_BUFFER,
	  PIXEL_UNPACK_BUFFER = GL_PIXEL_UNPACK_BUFFER
	};
    }

    /*! \brief The possible host access patterns, if the host
     * accesses the data.
     */
    namespace buffer_usage {
      enum Enum
	{
	  STREAM_DRAW = GL_STREAM_DRAW, 
	  STREAM_READ = GL_STREAM_READ, 
	  STREAM_COPY = GL_STREAM_COPY, 
	  STATIC_DRAW = GL_STATIC_DRAW, 
	  STATIC_READ = GL_STATIC_READ, 
	  STATIC_COPY = GL_STATIC_COPY, 
	  DYNAMIC_DRAW = GL_DYNAMIC_DRAW, 
	  DYNAMIC_READ = GL_DYNAMIC_READ, 
	  DYNAMIC_COPY = GL_DYNAMIC_COPY
	};
    }
  }
}
