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

namespace magnet {
  namespace GL {
    namespace objects {
      /*! \brief Base class for objects using instancing.
       *
       */
      class Instanced
      {	
      public:
	inline Instanced():_N(0) {}
	inline ~Instanced() { deinit(); }

	//! \brief Release any associated OpenGL resources.
	inline void deinit() 
	{
	  _N = 0;
	  _positionData.deinit();
	  _orientationData.deinit();
	  _colorData.deinit();
	  _primitiveVertices.deinit();
	  _primitiveNormals.deinit();
	  _primitiveIndices.deinit();
	}

	/*! \brief Initializes the OpenGL buffers.
	 *
	 * \param N The number of instances of the primitive object to draw
	 */
	inline void init(size_t N)
	{
	  deinit();
	  _N = N;

	  //Load the primitive data into the VBO's
	  _primitiveVertices.init(getPrimitiveVertices(), buffer_usage::STATIC_DRAW);
	  _primitiveNormals.init(getPrimitiveNormals(), buffer_usage::STATIC_DRAW);
	  _primitiveIndices.init(getPrimitiveIndicies(), buffer_usage::STATIC_DRAW);

	  {//Test positions!
	    std::vector<GLfloat> vertices(3 * _N);
	    
	    for (size_t i = 0; i < _N; ++i)
	      {
		vertices[3 * i + 0] = i * 2;
		vertices[3 * i + 1] = 0;
		vertices[3 * i + 2] = 0;
	      }
	    
	    _positionData.init(vertices);
	  }
	}
	
	/*! \brief Renders the instanced object.
	 */
	inline void glRender()
	{
	  _primitiveVertices.attachToVertex(3);
	  _primitiveNormals.attachToNormal();
	  glColor4f(0,1,1,1);
	  _primitiveIndices.drawInstancedElements(element_type::TRIANGLE_STRIP, _N);

	  glDisableVertexAttribArray(0);
	  glDisableClientState(GL_NORMAL_ARRAY);
	  glDisableClientState(GL_VERTEX_ARRAY);
	}

	virtual std::vector<GLfloat> getPrimitiveVertices() = 0;
	virtual std::vector<GLfloat> getPrimitiveNormals()  = 0;
	virtual std::vector<GLuint>  getPrimitiveIndicies() = 0;

      protected:
	size_t _N;
	magnet::GL::Buffer<GLfloat> _positionData;
	magnet::GL::Buffer<GLfloat> _orientationData;
	magnet::GL::Buffer<GLubyte> _colorData;
	magnet::GL::Buffer<GLfloat> _primitiveVertices;
	magnet::GL::Buffer<GLfloat> _primitiveNormals;
	magnet::GL::Buffer<GLuint>  _primitiveIndices;
      };
    }
  }
}
