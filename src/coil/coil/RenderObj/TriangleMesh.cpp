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

#include "TriangleMesh.hpp"
#include <magnet/math/vector.hpp>

namespace coil {
  void 
  RTriangleMesh::init(const magnet::thread::RefPtr<magnet::thread::TaskQueue>& systemQueue)
  {
    RTriangles::init(systemQueue);

    //Send the data we already have
    setGLPositions(_vertices);
    setGLElements(_elements);

    {//Calculate the normal vectors
      std::vector<float> VertexNormals(_vertices.size(), 0);
    
      //For every triangle, add the cross product of the two edges. We
      //then renormalize the normal to get a
      //"weighted-by-the-triangle-size" normal.
    
      for (size_t triangle(0); triangle < _elements.size() / 3; ++triangle)
	{
	  //Grab the vertex IDs
	  size_t v1(_elements[3 * triangle + 0]),
	    v2(_elements[3 * triangle + 1]),
	    v3(_elements[3 * triangle + 2]);
	
	  Vector V1(_vertices[3 * v1 + 0],
		    _vertices[3 * v1 + 1],
		    _vertices[3 * v1 + 2]),
	    V2(_vertices[3 * v2 + 0],
	       _vertices[3 * v2 + 1],
	       _vertices[3 * v2 + 2]),
	    V3(_vertices[3 * v3 + 0],
	       _vertices[3 * v3 + 1],
	       _vertices[3 * v3 + 2]);

	  Vector norm = (V2-V1)^(V3-V2);
	
	  for (size_t i(0); i < 3; ++i)
	    {
	      VertexNormals[3 * v1 + i] += norm[i];
	      VertexNormals[3 * v2 + i] += norm[i];
	      VertexNormals[3 * v3 + i] += norm[i];
	    }
	}

      //Now normalize those vertices
      for (size_t vert(0); vert < _vertices.size() / 3; ++vert)
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
  
    //Reclaim some memory
    _vertices.clear();
    _elements.clear();
  }
}
