/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

    This program is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    version 3 as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#pragma once
#include <magnet/GL/buffer.hpp>
#include <tr1/array>

namespace magnet {
  namespace GL {
    namespace objects {
      /*! \brief A regular grid object.
       *
       * This object is a regular grid with a set number of lines in
       * the x and y directions.
       *
       * This grid is centered on [0,0,0] and lies in
       * [\f$\pm0.5\f$,\f$\pm0.5\f$,0]. If you need the grid at
       * another location or with a different size then modify the
       * modelview matrix with scale and translate commands.
       */
      class Grid
      {
      public:
	//! \brief Destructor
	inline ~Grid() { deinit(); }

	//! \brief Release any associated OpenGL resources.
	inline void deinit() { _vertexData.deinit(); _xGridLines = _yGridLines = 0; }

	/*! \brief Sets up the vertex buffer objects for the regular
	 * grid.
	 *
	 * \param xlines The number of grid lines in the x dimension.
	 * \param ylines The number of grid lines in the y dimension.
	 */
	inline void init(size_t xlines, size_t ylines)
	{
	  _xGridLines = xlines;
	  _yGridLines = ylines;

	  std::vector<GLfloat> data(6 * (_xGridLines + _yGridLines));

	  for (size_t i(0); i < _xGridLines; ++i)
	    {
	      data[(i * 2 + 0) * 3 + 0] = -0.5f + i / float(_xGridLines - 1);
	      data[(i * 2 + 0) * 3 + 1] = -0.5f;
	      data[(i * 2 + 0) * 3 + 2] = 0;
	      data[(i * 2 + 1) * 3 + 0] = -0.5f + i / float(_xGridLines - 1);
	      data[(i * 2 + 1) * 3 + 1] = 0.5f;
	      data[(i * 2 + 1) * 3 + 2] = 0;
	    }
    
	  for (size_t i(0); i < _yGridLines; ++i)
	    {
	      data[((i + _xGridLines) * 2 + 0) * 3 + 0] = -0.5f;
	      data[((i + _xGridLines) * 2 + 0) * 3 + 1] = -0.5f + i / float(_yGridLines - 1);
	      data[((i + _xGridLines) * 2 + 0) * 3 + 2] = 0;
	      data[((i + _xGridLines) * 2 + 1) * 3 + 0] = 0.5f;
	      data[((i + _xGridLines) * 2 + 1) * 3 + 1] = -0.5f + i / float(_yGridLines - 1);
	      data[((i + _xGridLines) * 2 + 1) * 3 + 2] = 0;
	    }

	  _vertexData.init(data);
	}

	/*! \brief Attaches the vertex buffer and renders the regular grid.
	 *
	 * The color of the grid should be set before calling this
	 * function with glColor.
	 */
	inline void glRender()
	{
	  if (!(_xGridLines + _yGridLines))
	    M_throw() << "Cannot render uninitialized Grid object.";

	  _vertexData.bind(magnet::GL::ARRAY);
	  glVertexPointer(3, GL_FLOAT, 0, 0);
	  glEnableClientState(GL_VERTEX_ARRAY);

	  glDrawArrays(GL_LINES, 0, 6 * (_xGridLines + _xGridLines));

	  glDisableClientState(GL_VERTEX_ARRAY);
	}

      protected:
	magnet::GL::Buffer<GLfloat> _vertexData;
	size_t _xGridLines;
	size_t _yGridLines;
      };
    }
  }
}
