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
#include <magnet/GL/objects/instanced.hpp>

namespace magnet {
  namespace GL {
    namespace objects {
      /*! \brief A collection of cylinders.
       *
       */
      class Cylinders : public Instanced
      {	
      public:
	inline void init(size_t N, size_t LOD = 6)
	{
	  _LOD = LOD;
	  Instanced::init(N);
	}
	
	virtual element_type::Enum  getElementType() { return element_type::TRIANGLE_STRIP; }
	
	virtual std::vector<GLfloat> getPrimitiveVertices()
	{
	  std::vector<GLfloat> vertices(2 * _LOD * 3);

	  for (size_t vert = 0; vert < 2 * _LOD; ++vert)
	    {
	      vertices[3 * vert + 0] = 0.5f * std::sin((vert / 2) * 2.0f * M_PI / _LOD);
	      vertices[3 * vert + 1] = 0.5f * std::cos((vert / 2) * 2.0f * M_PI / _LOD);
	      vertices[3 * vert + 2] = vert % 2;
	    }
	  return vertices;
	}

	virtual std::vector<GLfloat> getPrimitiveNormals()
	{
	  std::vector<GLfloat> normals(2 * _LOD * 3);

	  for (size_t vert = 0; vert < 2 * _LOD; ++vert)
	    {
	      GLfloat x = 0.5f * std::sin((vert / 2) * 2.0f * M_PI / _LOD);
	      GLfloat y = 0.5f * std::cos((vert / 2) * 2.0f * M_PI / _LOD);
	      GLfloat scale = 1.0f / std::sqrt(x * x + y * y);
	      normals[3 * vert + 0] = x * scale;
	      normals[3 * vert + 1] = y * scale;
	      normals[3 * vert + 2] = 0;
	    }
	  return normals;
	}

	virtual std::vector<GLuint>  getPrimitiveIndicies()
	{
	  std::vector<GLuint> indices(2 * _LOD + 2);
	  
	  //Main vertices
	  for (size_t vert = 0; vert < 2 * _LOD; ++vert)
	    indices[vert] = vert;
	  
	  //Rejoin the end vertex
	  indices[2 * _LOD + 0] =  0;
	  indices[2 * _LOD + 1] =  1;
	  
	  return indices;
	}

      protected:
	size_t _LOD;
      };
    }
  }
}
