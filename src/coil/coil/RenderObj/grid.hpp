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

class Grid
{
public:
  inline ~Grid() { deinit(); }

  inline void deinit() { _renderData.deinit(); _xGridLines = _yGridLines = 0; }

  inline void init(size_t xlines, size_t ylines)
  {
    _xGridLines = xlines;
    _yGridLines = ylines;

    std::vector<GLfloat> data(6 * (_xGridLines + _yGridLines + 2));

    for (size_t i(0); i <= _xGridLines; ++i)
      {
	data[(i * 2 + 0) * 3 + 0] = -0.5f + i / float(numGridlines);
	data[(i * 2 + 0) * 3 + 1] = -0.5f;
	data[(i * 2 + 0) * 3 + 2] = 0;
	data[(i * 2 + 1) * 3 + 0] = -0.5f + i / float(numGridlines);
	data[(i * 2 + 1) * 3 + 1] = 0.5f;
	data[(i * 2 + 1) * 3 + 2] = 0;
      }
    
    for (size_t i(_xGridLines + 1); i <= _xGridLines + _yGridLines + 1; ++i)
      {
	data[(i * 2 + 0) * 3 + 0] = -0.5f;
	data[(i * 2 + 0) * 3 + 1] = -0.5f + i / float(numGridlines);
	data[(i * 2 + 0) * 3 + 2] = 0;
	data[(i * 2 + 1) * 3 + 0] = 0.5f;
	data[(i * 2 + 1) * 3 + 1] = -0.5f + i / float(numGridlines);
	data[(i * 2 + 1) * 3 + 2] = 0;
      }

    _renderData.init(data);
  }

  inline void glRender()
  {
    if (!(_xGridLines + _yGridLines))
      M_throw() << "Cannot render uninitialized Grid object.";

    _renderData.bind(magnet::GL::Buffer::ARRAY);
    glVertexPointer(3, GL_FLOAT, 0, 0);
    glEnableClientState(GL_VERTEX_ARRAY);
    glDrawArrays(GL_LINES, 0, 6 * (_xGridLines + _xGridLines + 2));
    glDisableClientState(GL_VERTEX_ARRAY);
  }

protected:
  magnet::GL::Buffer _renderData;
  size_t _xGridLines;
  size_t _yGridLines;
};
