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

#include <coil/RenderObj/Surface.hpp>
#include <magnet/string/formatcode.hpp>
#include <magnet/string/line_number.hpp>
#include <coil/images/images.hpp>
#include <iostream>

namespace coil {
  Glib::RefPtr<Gdk::Pixbuf> 
  RSurface::getIcon()
  { return coil::images::Function_Icon(); }

  RSurface::RSurface(std::string name,
		     size_t N, Vector origin, Vector axis1,
		     Vector axis2, Vector axis3):
    RTriangles(name),
    _N(N),
    _origin(origin),
    _axis1(axis1),
    _axis2(axis2),
    _axis3(axis3)
  {
    //Ensure that the axis have the correct ordering (so the normals
    //aren't messed up as we use CW/CCW vertex ordering).
    if (((axis1 ^ axis2) | axis3) > 0) std::swap(_axis1, _axis2);
  }

  void 
  RSurface::init(const std::shared_ptr<magnet::thread::TaskQueue>& systemQueue)
  {
    RTriangles::init(systemQueue);

    {//Setup initial vertex positions
      std::vector<float> VertexPos(3 * _N * _N, 0.0);
      
      Vector axis1step = _axis1 / (_N-1), axis2step = _axis2 / (_N-1);

      for (size_t i = 0; i < _N; ++i)
	for (size_t j = 0; j < _N; ++j)
	  {
	    Vector pos = i * axis1step + j * axis2step + _origin;
	    VertexPos[0 + 3 * (i + _N * j)] = pos[0];	
	    VertexPos[1 + 3 * (i + _N * j)] = pos[1];
	    VertexPos[2 + 3 * (i + _N * j)] = pos[2];
	  }

      setGLPositions(VertexPos);
    }

    {//Setup inital normal vectors
      std::vector<float> VertexNormals(3 * _N * _N, 0.0);
      Vector normal = _axis3 / _axis3.nrm();

      for (size_t i = 0; i < _N; ++i)
	for (size_t j = 0; j < _N; ++j)
	  {
	    VertexNormals[0 + 3 * (i + _N * j)] = normal[0];
	    VertexNormals[1 + 3 * (i + _N * j)] = normal[1];
	    VertexNormals[2 + 3 * (i + _N * j)] = normal[2];
	  }
      setGLNormals(VertexNormals);
    }

    {//Setup initial Colors
      std::vector<GLubyte> VertexColor(_N * _N * 4);
      for (size_t i = 0; i < _N; ++i)
	for (size_t j = 0; j < _N; ++j)
	  {
	    VertexColor[(i + _N * j) * 4 + 0] = 255;
	    VertexColor[(i + _N * j) * 4 + 1] = 255;
	    VertexColor[(i + _N * j) * 4 + 2] = 255;
	    VertexColor[(i + _N * j) * 4 + 3] = 255;
	  }
      setGLColors(VertexColor);
    }
   
    {//Setup initial element data
      std::vector<GLuint> ElementData(3*2*(_N-1)*(_N-1), 0.0);
      for  (size_t i = 0; i < _N - 1; i++)
	for (size_t j = 0; j < _N - 1; j++)
	  {
	    ElementData[6 * (i + (_N - 1) * j) + 0] = i + _N * j;
	    ElementData[6 * (i + (_N - 1) * j) + 1] = i + _N * (j + 1);
	    ElementData[6 * (i + (_N - 1) * j) + 2] = i + 1 + _N * (j + 1);
	    ElementData[6 * (i + (_N - 1) * j) + 3] = i + _N * j;
	    ElementData[6 * (i + (_N - 1) * j) + 4] = i + 1 + _N * (j + 1);
	    ElementData[6 * (i + (_N - 1) * j) + 5] = i + 1 + _N * j;
	  }
      setGLElements(ElementData);
    }

    _initialised = True;
  }

  magnet::math::Vector 
  RSurface::getMaxCoord() const
  {
    magnet::math::Vector max = _origin;
    for (size_t i(0); i< 3; ++i)
      {
	max[i] = std::max(_origin[i] + _axis1[i], max[i]);
	max[i] = std::max(_origin[i] + _axis2[i], max[i]);
      }

    return max;
  }

  magnet::math::Vector 
  RSurface::getMinCoord() const
  {
    magnet::math::Vector min = _origin;
    for (size_t i(0); i< 3; ++i)
      {
	min[i] = std::min(_origin[i] + _axis1[i], min[i]);
	min[i] = std::min(_origin[i] + _axis2[i], min[i]);
      }

    return min;
  }
}
