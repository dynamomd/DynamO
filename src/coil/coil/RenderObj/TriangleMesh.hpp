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

#include "Triangles.hpp"

namespace coil {
  class RTriangleMesh : public RTriangles
  {
  public:
    RTriangleMesh(std::string name, 
		  std::vector<GLfloat>& vertices,
		  std::vector<GLuint>& elements):
      RTriangles(name),
      _vertices(vertices),
      _elements(elements)
    {}

    virtual void init(const magnet::thread::RefPtr<magnet::thread::TaskQueue>& systemQueue);

  protected:
    std::vector<GLfloat> _vertices;
    std::vector<GLuint> _elements;
  };
}

