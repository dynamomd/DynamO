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
#include <vector>
#include <cmath>

#include <string.h>
#include <stdlib.h>
#include <magnet/math/vector.hpp>
#include <magnet/exception.hpp>

namespace magnet {
  namespace GL {
    namespace objects {
      namespace primitives {
	/*! \brief This class contains functions which generate the
            vertex data for an OpenGL cylinder.
	 */
	class Cylinder
	{
	public:

	  inline static std::vector<GLfloat> getVertices(size_t LOD)
	  {
	    std::vector<GLfloat> vertices(2 * LOD * 3);
	    
	    for (size_t vert = 0; vert < 2 * LOD; ++vert)
	      {
		vertices[3 * vert + 0] = 0.5f * std::sin((vert / 2) * 2.0f * M_PI / LOD);
		vertices[3 * vert + 1] = 0.5f * std::cos((vert / 2) * 2.0f * M_PI / LOD);
		vertices[3 * vert + 2] = vert % 2;
	      }
	    return vertices;
	  }

	  inline static std::vector<GLfloat> getNormals(size_t LOD)
	  {
	    std::vector<GLfloat> normals(2 * LOD * 3);

	    for (size_t vert = 0; vert < 2 * LOD; ++vert)
	      {
		GLfloat x = 0.5f * std::sin((vert / 2) * 2.0f * M_PI / LOD);
		GLfloat y = 0.5f * std::cos((vert / 2) * 2.0f * M_PI / LOD);
		GLfloat scale = 1.0f / std::sqrt(x * x + y * y);
		normals[3 * vert + 0] = x * scale;
		normals[3 * vert + 1] = y * scale;
		normals[3 * vert + 2] = 0;
	      }
	    return normals;
	  }

	  inline static std::vector<GLuint> getIndices(size_t LOD)
	  {
	    //2 triangles per face (3 indices per triangle)
	    std::vector<GLuint> indices(6 * LOD);      
	    for (size_t vert = 0; vert < LOD; ++vert)
	      {
		indices[6 * vert + 0] = (2 * vert + 0) % (2 * LOD);
		indices[6 * vert + 1] = (2 * vert + 1) % (2 * LOD);
		indices[6 * vert + 2] = (2 * vert + 2) % (2 * LOD);
		indices[6 * vert + 3] = (2 * vert + 1) % (2 * LOD);
		indices[6 * vert + 4] = (2 * vert + 3) % (2 * LOD);
		indices[6 * vert + 5] = (2 * vert + 2) % (2 * LOD);
	      }
	    return indices;
	  }

	protected:
	};
      }
    }
  }
}
