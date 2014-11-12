/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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

#include "TriangleMesh.hpp"
#include <magnet/math/vector.hpp>

namespace coil {
  void 
  RTriangleMesh::init(const std::shared_ptr<magnet::thread::TaskQueue>& systemQueue)
  {
    RTriangles::init(systemQueue);
    //_context = magnet::GL::Context::getContext();
    updateGLDataWorker(_vertices, _elements, _colours);
    _vertices.clear();
    _elements.clear();
    _colours.clear();
  }

  void 
  RTriangleMesh::updateGLDataWorker(const std::vector<GLfloat> vertices, const std::vector<GLuint> elements, const std::vector<GLubyte> colours)
  {
    //Send the data we already have
    setGLPositions(vertices);
    setGLElements(elements);
    {//Calculate the normal vectors                          
      std::vector<float> VertexNormals(3 * vertices.size() / _triangleComponents, 0);
    
      //For every triangle, add the cross product of the two edges. We
      //then renormalize the normal to get a
      //"weighted-by-the-triangle-size" normal.
    
      for (size_t triangle(0); triangle < elements.size() / 3; ++triangle)
	{
	  //Grab the vertex IDs
	  size_t v1(elements[3 * triangle + 0]),
	    v2(elements[3 * triangle + 1]),
	    v3(elements[3 * triangle + 2]);
	  
	  Vector V1{vertices[3 * v1 + 0], vertices[3 * v1 + 1], vertices[3 * v1 + 2]},
	    V2{vertices[3 * v2 + 0], vertices[3 * v2 + 1], vertices[3 * v2 + 2]},
	      V3{vertices[3 * v3 + 0], vertices[3 * v3 + 1], vertices[3 * v3 + 2]};

	  Vector norm = (V2-V1)^(V3-V2);
	
	  for (size_t i(0); i < 3; ++i)
	    {
	      VertexNormals[3 * v1 + i] += norm[i];
	      VertexNormals[3 * v2 + i] += norm[i];
	      VertexNormals[3 * v3 + i] += norm[i];
	    }
	}

      //Now normalize those vertices
      for (size_t vert(0); vert < vertices.size() / 3; ++vert)
	{
	  double norm 
	    = VertexNormals[3 * vert + 0] * VertexNormals[3 * vert + 0]
	    + VertexNormals[3 * vert + 1] * VertexNormals[3 * vert + 1]
	    + VertexNormals[3 * vert + 2] * VertexNormals[3 * vert + 2];
	
	  if (norm)
	    {
	      double factor = 1 / std::sqrt(norm);
	      VertexNormals[3 * vert + 0] *= factor;
	      VertexNormals[3 * vert + 1] *= factor;
	      VertexNormals[3 * vert + 2] *= factor;
	    }
	  else
	    {
	      VertexNormals[3 * vert + 0] = 1;
	      VertexNormals[3 * vert + 1] = 0;
	      VertexNormals[3 * vert + 2] = 0;
	    }
	}

      setGLNormals(VertexNormals);
    }

    if (colours.empty())
      setGLColors(std::vector<GLubyte>((vertices.size() / _triangleComponents) * 4, 255));
    else
      setGLColors(colours);
  }
}
