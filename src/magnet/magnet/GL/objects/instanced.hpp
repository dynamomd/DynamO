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
	  //No need to deinitialise, we'll just initialise over the top of the old data
	  //deinit(); 
	  _N = N;
	  
	  //Load the primitive data into the VBO's
	  _primitiveVertices.init(getPrimitiveVertices(), buffer_usage::STATIC_DRAW);
	  _primitiveNormals.init(getPrimitiveNormals(), buffer_usage::STATIC_DRAW);
	  _primitiveIndices.init(getPrimitiveIndicies(), buffer_usage::STATIC_DRAW);
	}

	/*! \brief Renders the instanced object.
	 */
	inline void glRender()
	{
	  _primitiveVertices.attachToVertex();
	  _primitiveNormals.attachToNormal();
	  _primitiveIndices.drawInstancedElements(getElementType(), _N);
	}

	virtual std::vector<GLfloat> getPrimitiveVertices() = 0;
	virtual std::vector<GLfloat> getPrimitiveNormals()  = 0;
	virtual std::vector<GLuint>  getPrimitiveIndicies() = 0;
	virtual element_type::Enum   getElementType() = 0;

      protected:
	size_t _N;
	magnet::GL::Buffer<GLfloat> _primitiveVertices;
	magnet::GL::Buffer<GLfloat> _primitiveNormals;
	magnet::GL::Buffer<GLuint>  _primitiveIndices;
      };
    }
  }
}
