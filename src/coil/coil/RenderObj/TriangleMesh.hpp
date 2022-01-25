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

#include "Triangles.hpp"

namespace coil {
  class RTriangleMesh : public RTriangles
  {
  public:
    RTriangleMesh(magnet::GL::Context::ContextPtr context, std::string name, 
		  const std::vector<GLfloat>& vertices,
		  const std::vector<GLuint>& elements,
		  const std::vector<GLubyte>& colours = std::vector<GLubyte>()):
      RTriangles(context, name),
      _vertices(vertices),
      _elements(elements),
      _colours(colours)
    {}

    virtual void init(const std::shared_ptr<magnet::thread::TaskQueue>& systemQueue);

    void updateGLData(const std::vector<GLfloat> vertices, const std::vector<GLuint> elements, const std::vector<GLubyte> colours = std::vector<GLubyte>()) {
      _context->queueTask(std::bind(&RTriangleMesh::updateGLDataWorker, this, vertices, elements, colours));    
    }

  protected:
    void updateGLDataWorker(const std::vector<GLfloat> vertices, const std::vector<GLuint> elements, const std::vector<GLubyte> colours = std::vector<GLubyte>());

    std::vector<GLfloat> _vertices;
    std::vector<GLuint> _elements;
    std::vector<GLubyte> _colours;    
    magnet::GL::Context::ContextPtr _context;
  };
}

