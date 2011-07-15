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
#include "RenderObj.hpp"
#include <magnet/GL/buffer.hpp>
#include <vector>


class RLines : public RenderObj
{
public:
  RLines(size_t N, std::string name);
  ~RLines();

  virtual void glRender();
  virtual void initOpenGL();

  void setGLColors(std::vector<GLubyte>& VertexColor);
  void setGLPositions(std::vector<GLfloat>& VertexPos);
  void setGLElements(std::vector<GLuint>& Elements);

  magnet::GL::Buffer<GLfloat>& getVertexGLData() { return _posBuff; }
  magnet::GL::Buffer<GLubyte>& getColorGLData() { return _colBuff; }

  virtual void releaseCLGLResources();

protected:
  size_t _N;
  magnet::GL::Buffer<GLubyte> _colBuff;
  magnet::GL::Buffer<GLfloat> _posBuff;  
  magnet::GL::Buffer<GLuint> _elementBuff;
};
