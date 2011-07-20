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
      /*! \brief A simple cube object.
       */
      class Cube
      {
      public:
	//! \brief Destructor
	inline ~Cube() { deinit(); }

	//! \brief Release any associated OpenGL resources.
	inline void deinit() { _vertexData.deinit(); _indexData.deinit(); }

	/*! \brief Sets up the vertex buffer objects for the cube object.
	 */
	inline void init()
	{
	  GLfloat vertices[] = {-1,-1,-1,  1,-1,-1,  1, 1,-1, -1, 1,-1,
				-1,-1, 1, -1, 1, 1,  1, 1, 1,  1,-1, 1};
	  GLuint elements[] = {3,2,1,0, 6,7,1,2, 5,4,7,6, 3,0,4,5, 6,2,3,5, 7,4,0,1 };
	  
	  _vertexData.init(std::vector<GLfloat>(vertices, vertices + sizeof(vertices) / sizeof(float)));
	  _indexData.init(std::vector<GLuint>(elements, elements + sizeof(elements) / sizeof(int)));
	}

	/*! \brief Attaches the vertex buffer and renders the regular
	 * cube.
	 */
	inline void glRender()
	{
	  _vertexData.attachToVertex();
	  _indexData.drawElements(magnet::GL::element_type::QUADS);
	}

      protected:
	magnet::GL::Buffer<GLfloat> _vertexData;
	magnet::GL::Buffer<GLuint>  _indexData;
      };
    }
  }
}
