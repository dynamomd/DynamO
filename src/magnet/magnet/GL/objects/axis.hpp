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
#include <tr1/array>

namespace magnet {
  namespace GL {
    namespace objects {
      /*! \brief An axis for indicating the orientation of a render.
       *
       * The axis is centered on [0,0,0] and lies in
       * [\f$\pm0.5\f$,\f$\pm0.5\f$,\f$\pm0.5\f$]. If you need the axis at
       * another location or with a different size then modify the
       * modelview matrix with scale and translate commands.
       */
      class Axis
      {
      public:
	//! \brief Destructor
	inline ~Axis() { deinit(); }

	//! \brief Release any associated OpenGL resources.
	inline void deinit() { _vertexData.deinit(); _colorData.deinit(); }

	/*! \brief Sets up the vertex buffer objects for the axis.
	 */
	inline void init()
	{
	  ///////////////////Vertex Data
	  // three arrows, consisting of 3 lines each (6 vertices
	  std::vector<GLfloat> vertexdata(18 * 3);

	  float arrowHeadDepth = 0.25f;
	  float arrowHeadWidth = 0.15f;

	  for (size_t arrow(0); arrow < 3; ++arrow)
	    {
	      //Ordinate ->
	      for (size_t i(0); i < 3; ++i)
		vertexdata[(6 * arrow + 0) * 3 + i] = -0.5f;

	      //-> Head
	      for (size_t i(0); i < 3; ++i)
		vertexdata[(6 * arrow + 1) * 3 + i] = -0.5f;
	      vertexdata[(6 * arrow + 1) * 3 + arrow] = 0.5f;

	      //Head ->
	      for (size_t i(0); i < 3; ++i)
		vertexdata[(6 * arrow + 2) * 3 + i] = -0.5f;
	      vertexdata[(6 * arrow + 2) * 3 + arrow] = 0.5f;

	      //-> Left point
	      vertexdata[(6 * arrow + 3) * 3 + arrow] = 0.5f - arrowHeadDepth;
	      vertexdata[(6 * arrow + 3) * 3 + ((arrow + 1) % 3)] = -0.5f + arrowHeadWidth;
	      vertexdata[(6 * arrow + 3) * 3 + ((arrow + 2) % 3)] = -0.5f;

	      //Head ->
	      for (size_t i(0); i < 3; ++i)
		vertexdata[(6 * arrow + 4) * 3 + i] = -0.5f;
	      vertexdata[(6 * arrow + 4) * 3 + arrow] = 0.5f;

	      //-> Right point
	      vertexdata[(6 * arrow + 5) * 3 + arrow] = 0.5f - arrowHeadDepth;
	      vertexdata[(6 * arrow + 5) * 3 + ((arrow + 1) % 3)] = -0.5f - arrowHeadWidth;
	      vertexdata[(6 * arrow + 5) * 3 + ((arrow + 2) % 3)] = -0.5f;
	    }

	  _vertexData.init(vertexdata);
	  ///////////////////Color Data
	  std::vector<GLubyte> colordata(18 * 4);
	  
	  //First arrow is red, second green and last blue
	  for (size_t arrow(0); arrow < 3; ++arrow)
	    for (size_t vertex(0); vertex < 6; ++vertex)
	      {
		colordata[(arrow * 6 + vertex)* 4 + arrow] = 255;
		colordata[(arrow * 6 + vertex)* 4 + ((arrow + 1) % 3)] = 0;
		colordata[(arrow * 6 + vertex)* 4 + ((arrow + 2) % 3)] = 0;
		colordata[(arrow * 6 + vertex)* 4 + 3] = 255; //Alpha channel
	      }

	  _colorData.init(colordata);
	}

	/*! \brief Attaches the vertex buffer and renders the axis.
	 */
	inline void glRender()
	{
	  _colorData.getContext()->cleanupAttributeArrays();
	  _colorData.attachToColor();
	  _vertexData.drawArray(magnet::GL::element_type::LINES);
	}

	const Context::ContextPtr& getContext() { return _vertexData.getContext(); }
      protected:
	magnet::GL::Buffer<GLfloat> _vertexData;
	magnet::GL::Buffer<GLubyte> _colorData;
      };
    }
  }
}
