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
      /*! \brief A collection of cylinders.
       *
       */
      class Cylinders
      {
      public:
	//! \brief Destructor
	inline ~Cylinders() { deinit(); }

	//! \brief Release any associated OpenGL resources.
	inline void deinit() { _vertexData.deinit(); _indexData.deinit(); _cylinders = _LOD = 0; }

	/*! \brief Sets up the vertex buffer objects for the regular
	 * grid.
	 *
	 * \param LOD The number of vertices on the edge of the
	 * cylinder (Level of Detail).
	 *
	 * \param cylinders The number of cylinders represented by this object
	 */
	inline void init(size_t cylinders, size_t LOD = 6)
	{
	  deinit();
	  _LOD = LOD;
	  _cylinders = cylinders;

	  //Every cylinder is made of 2 * LOD vertices, which loop back
	  //on itself (+2) and have a degenerate index before and
	  //after (+2)

	  {
	    std::vector<GLuint> indices((2 * LOD + 4) * _cylinders);
	    
	    for (size_t cyl(0); cyl < _cylinders; ++cyl)
	      {
		//Degenerate "begin" vertex
		indices[(2 * LOD + 4) * cyl + 0] =  (2 * LOD) * _cylinders + 0;
		
		//Main vertices
		for (size_t vert = 0; vert < 2 * LOD; ++vert)
		  indices[(2 * LOD + 4) * cyl + 1 + vert] =  (2 * LOD) * _cylinders + vert;
		
		//Rejoin the end vertex
		indices[(2 * LOD + 4) * cyl + 2 * LOD + 1] =  (2 * LOD) * _cylinders + 0;
		indices[(2 * LOD + 4) * cyl + 2 * LOD + 2] =  (2 * LOD) * _cylinders + 1;
		
		//Degenerate "end" vertex
		indices[(2 * LOD + 4) * cyl + 2 * LOD + 3] =  (2 * LOD) * _cylinders + 1;
	      }
	    
	    _indexData.init(indices, GL_STATIC_DRAW);
	  }
	  
	  //Reserve the vertex data buffer
	  _vertexData.init(2 * LOD * _cylinders * 3 * sizeof(GLfloat));
	}

	/*! \brief Attaches the vertex buffer and renders the regular grid.
	 *
	 * The color of the grid should be set before calling this
	 * function with glColor.
	 */
	inline void glRender()
	{
	  if (!_cylinders)
	    M_throw() << "Trying to render an empty set of cylinders";

	  _vertexData.bind(magnet::GL::Buffer::ARRAY);
	  glVertexPointer(3, GL_FLOAT, 0, 0);
	  glEnableClientState(GL_VERTEX_ARRAY);

	  _elementBuff.bind(magnet::GL::Buffer::ELEMENT_ARRAY);
	  glDrawElements(GL_TRIANGLE_STRIP, (2 * LOD + 4) * _cylinders, GL_UNSIGNED_INT, 0);

	  glDisableClientState(GL_VERTEX_ARRAY);
	}

      protected:
	magnet::GL::Buffer<GLfloat> _vertexData;
	magnet::GL::Buffer<GLuint> _indexData;
	size_t _LOD;
	size_t _cylinders;
      };
    }
  }
}
