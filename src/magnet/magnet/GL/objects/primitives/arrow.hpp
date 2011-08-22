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
	class Arrow
	{
	public:

	  inline static std::vector<GLfloat> getVertices(size_t LOD, GLfloat head_length_ratio = 0.5,
							 GLfloat body_radius_ratio = 0.5)
	  {
	    //4 body vertices per LOD
	    std::vector<GLfloat> vertices;
	    for (size_t slice = 0; slice < LOD; ++slice)
	      {		
		GLfloat angle = slice * 2.0f * M_PI / LOD;
		GLfloat x = 0.25f * std::sin(angle);
		GLfloat y = 0.25f * std::cos(angle);

		//Add the point vertex for the arrow
		vertices.push_back(0); 
		vertices.push_back(0); 
		vertices.push_back(1);
		
		//Head cone vertex
		vertices.push_back(x);
		vertices.push_back(y);
		vertices.push_back(1 - head_length_ratio);

		//Head-cylinder vertex
		vertices.push_back(x * body_radius_ratio);
		vertices.push_back(y * body_radius_ratio);
		vertices.push_back(1 - head_length_ratio);

		//tail-cylinder vertex
		vertices.push_back(x * body_radius_ratio);
		vertices.push_back(y * body_radius_ratio);
		vertices.push_back(0);

		//end-cylinder vertex
		vertices.push_back(x * body_radius_ratio);
		vertices.push_back(y * body_radius_ratio);
		vertices.push_back(0);

		//end-cylinder center vertex
		vertices.push_back(0);
		vertices.push_back(0);
		vertices.push_back(0);

		//Cone-cylinder outer vertex
		vertices.push_back(x);
		vertices.push_back(y);
		vertices.push_back(1 - head_length_ratio);

		//Cone-cylinder inner vertex
		vertices.push_back(x * body_radius_ratio);
		vertices.push_back(y * body_radius_ratio);
		vertices.push_back(1 - head_length_ratio);

	      }

	    return vertices;
	  }

	  inline static std::vector<GLfloat> getNormals(size_t LOD)
	  {
	    std::vector<GLfloat> normals;

	    for (size_t slice = 0; slice < LOD; ++slice)
	      {		
		GLfloat angle = slice * 2.0f * M_PI / LOD;
		GLfloat x = std::sin(angle);
		GLfloat y = std::cos(angle);

		math::Vector zaxis(0,0,1), radialaxis(x,y,0), edge(x,y,-1);

		edge /= edge.nrm();

		math::Vector rotationAxis = radialaxis ^ zaxis;

		math::Vector coneNormal = Rodrigues(rotationAxis * (M_PI / 2)) * edge;
		
		//Point vertex
		normals.push_back(coneNormal[0]); 
		normals.push_back(coneNormal[1]); 
		normals.push_back(coneNormal[2]);

		//Head cone vertex
		normals.push_back(coneNormal[0]); 
		normals.push_back(coneNormal[1]); 
		normals.push_back(coneNormal[2]);

		//Head-cylinder vertex
		normals.push_back(x);
		normals.push_back(y);
		normals.push_back(0);

		//tail-cylinder vertex
		normals.push_back(x);
		normals.push_back(y);
		normals.push_back(0);

		//end-cylinder vertex
		normals.push_back(0);
		normals.push_back(0);
		normals.push_back(-1);

		//end-cylinder center vertex
		normals.push_back(0);
		normals.push_back(0);
		normals.push_back(-1);

		//Cone-cylinder outer vertex
		normals.push_back(0);
		normals.push_back(0);
		normals.push_back(-1);

		//Cone-cylinder inner vertex
		normals.push_back(0);
		normals.push_back(0);
		normals.push_back(-1);

	      }

	    return normals;
	  }

	  inline static std::vector<GLuint> getIndices(size_t LOD)
	  {
	    //The number of vertices per body
	    size_t v_body = 8;
	    
	    //This is the number of body vertices
	    size_t vertex_count = LOD * v_body;

	    std::vector<GLuint> indices;

	    for (size_t slice = 0; slice < LOD; ++slice)
	      {
		//The cone triangle is first
		indices.push_back((v_body * (slice + 0) + 0) % vertex_count);
		indices.push_back((v_body * (slice + 1) + 1) % vertex_count);
		indices.push_back((v_body * (slice + 0) + 1) % vertex_count);

		//The body triangles
		indices.push_back((v_body * (slice + 1) + 2) % vertex_count);
		indices.push_back((v_body * (slice + 0) + 3) % vertex_count);
		indices.push_back((v_body * (slice + 0) + 2) % vertex_count);

		indices.push_back((v_body * (slice + 1) + 2) % vertex_count);
		indices.push_back((v_body * (slice + 1) + 3) % vertex_count);
		indices.push_back((v_body * (slice + 0) + 3) % vertex_count);

		//The cylinder end cap triangle
		indices.push_back((v_body * (slice + 0) + 4) % vertex_count);
		indices.push_back((v_body * (slice + 1) + 4) % vertex_count);
		indices.push_back((v_body * (slice + 0) + 5) % vertex_count);
		
		//The cone-cylinder connection triangles
		indices.push_back((v_body * (slice + 0) + 7) % vertex_count);
		indices.push_back((v_body * (slice + 0) + 6) % vertex_count);
		indices.push_back((v_body * (slice + 1) + 6) % vertex_count);

		indices.push_back((v_body * (slice + 0) + 7) % vertex_count);
		indices.push_back((v_body * (slice + 1) + 6) % vertex_count);
		indices.push_back((v_body * (slice + 1) + 7) % vertex_count);
	      }

	    return indices;
	  }

	protected:
	};
      }
    }
  }
}
