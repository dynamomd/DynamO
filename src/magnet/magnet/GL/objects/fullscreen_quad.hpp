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
      /*! \brief A full-screen quad.
       *
       * This object is used to generate a fragment shader for every
       * pixel on the screen. It assumes that the fragment shader is a
       * simple passthrough shader as the verticies are already in
       * screen/eye space.
       */
      class FullScreenQuad
      {
      public:
	//! \brief Destructor
	inline ~FullScreenQuad() { deinit(); }

	//! \brief Release any associated OpenGL resources.
	inline void deinit() { _vertexData.deinit(); }

	/*! \brief Sets up the vertex buffer objects for the quad.
	 */
	inline void init()
	{
	  ///////////////////Vertex Data
	  // Single quad, in pre-transformed screen coordinates
	  std::vector<GLfloat> vertexdata(4 * 2);
	  vertexdata[2 * 0 + 0] = -1; vertexdata[2 * 0 + 1] = -1;
	  vertexdata[2 * 1 + 0] =  1; vertexdata[2 * 1 + 1] = -1;
	  vertexdata[2 * 2 + 0] =  1; vertexdata[2 * 2 + 1] =  1;
	  vertexdata[2 * 3 + 0] = -1; vertexdata[2 * 3 + 1] =  1;
	  _vertexData.init(vertexdata);
	}

	/*! \brief Attaches the vertex buffer and renders the quad.
	 */
	inline void glRender()
	{ _vertexData.drawArray(magnet::GL::element_type::QUADS, 2); }

      protected:
	magnet::GL::Buffer<GLfloat> _vertexData;
      };
    }
  }
}
