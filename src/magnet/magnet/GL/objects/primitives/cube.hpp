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
#pragma once
#include <vector>

namespace magnet {
  namespace GL {
    namespace objects {
      namespace primitives {
	/*! \brief This class contains functions which generate the
            vertex data for an OpenGL cube.

	    An alternative render setup for a cube is the following,
	    but its use will generate incorrect face normals.

	  GLfloat vertices[] = {-0.5,-0.5,-0.5,  0.5,-0.5,-0.5,  0.5, 0.5,-0.5, -0.5, 0.5,-0.5,
				-0.5,-0.5, 0.5, -0.5, 0.5, 0.5,  0.5, 0.5, 0.5,  0.5,-0.5, 0.5};
	  GLuint elements[] = {3,2,1, 1,0,3, 
			       6,7,1, 1,2,6, 
			       5,4,7, 7,6,5, 
			       3,0,4, 4,5,3, 
			       6,2,3, 3,5,6, 
			       7,4,0, 0,1,7 };

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
