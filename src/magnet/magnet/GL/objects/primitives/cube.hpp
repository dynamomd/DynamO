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
            vertex data for an OpenGL cube.
	 */
	class Cube
	{
	public:

	  inline static std::vector<GLfloat> getVertices()
	  {
	    GLfloat vertices[] = {-0.5,0.5,-0.5, 0.5,0.5,-0.5, 0.5,-0.5,-0.5,
				  0.5,-0.5,-0.5, -0.5,-0.5,-0.5, -0.5, 0.5,-0.5,
				  0.5,0.5,0.5, 0.5,-0.5, 0.5, 0.5,-0.5,-0.5,
				  0.5,-0.5,-0.5, 0.5, 0.5,-0.5, 0.5,0.5,0.5,
				  -0.5,0.5,0.5, -0.5,-0.5,0.5, 0.5,-0.5,0.5,
				  0.5,-0.5,0.5, 0.5,0.5,0.5, -0.5,0.5,0.5,
				  -0.5, 0.5,-0.5, -0.5,-0.5,-0.5, -0.5,-0.5, 0.5,
				  -0.5,-0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5,-0.5,
				  0.5, 0.5, 0.5, 0.5, 0.5,-0.5, -0.5, 0.5,-0.5,
				  -0.5, 0.5,-0.5, -0.5, 0.5, 0.5,  0.5, 0.5, 0.5,
				  0.5,-0.5, 0.5,  -0.5,-0.5, 0.5, -0.5,-0.5,-0.5,
				  -0.5,-0.5,-0.5,  0.5,-0.5,-0.5, 0.5,-0.5, 0.5};

	    return std::vector<GLfloat>(vertices, vertices + sizeof(vertices) / sizeof(GLfloat));
	  }

	  inline static std::vector<GLfloat> getNormals()
	  {
	    GLfloat normals[] = {0,0,-1, 0,0,-1, 0,0,-1,
				 0,0,-1, 0,0,-1, 0,0,-1,
				 1,0,0, 1,0,0, 1,0,0,
				 1,0,0, 1,0,0, 1,0,0,
				 0,0,1, 0,0,1, 0,0,1,
				 0,0,1, 0,0,1, 0,0,1,
				 -1,0,0, -1,0,0, -1,0,0,
				 -1,0,0, -1,0,0, -1,0,0,
				 0,1,0, 0,1,0, 0,1,0,
				 0,1,0, 0,1,0, 0,1,0,
				 0,-1,0, 0,-1,0, 0,-1,0,
				 0,-1,0, 0,-1,0, 0,-1,0};
	    return std::vector<GLfloat>(normals, normals + sizeof(normals) / sizeof(GLfloat));
	  }

	  inline static std::vector<GLuint> getIndices()
	  {
	    std::vector<GLuint> retval; retval.resize(36);
	    for (size_t i(0); i < 36; ++i) retval[i]=i;
	    return retval;
	  }
	};
      }
    }
  }
}
