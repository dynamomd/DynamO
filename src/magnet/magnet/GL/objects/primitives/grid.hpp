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
/*! \brief This class generates the vertices for the lines of a regular grid.

  This grid is centered on [0,0,0] and lies in
  [\f$\pm0.5\f$,\f$\pm0.5\f$,0]. If you need the grid at
  another location or with a different size then modify the
  modelview matrix with scale and translate commands or edit
  the vertices by hand.
 */
class Grid {
public:
  inline static std::vector<GLfloat> getVertices(const size_t xlines,
                                                 const size_t ylines) {
    std::vector<GLfloat> data(6 * (xlines + ylines));

    for (size_t i(0); i < xlines; ++i) {
      data[(i * 2 + 0) * 3 + 0] = -0.5f + i / GLfloat(xlines - 1);
      data[(i * 2 + 0) * 3 + 1] = -0.5f;
      data[(i * 2 + 0) * 3 + 2] = 0;
      data[(i * 2 + 1) * 3 + 0] = -0.5f + i / GLfloat(ylines - 1);
      data[(i * 2 + 1) * 3 + 1] = 0.5f;
      data[(i * 2 + 1) * 3 + 2] = 0;
    }

    for (size_t i(0); i < ylines; ++i) {
      data[((i + xlines) * 2 + 0) * 3 + 0] = -0.5f;
      data[((i + xlines) * 2 + 0) * 3 + 1] = -0.5f + i / float(ylines - 1);
      data[((i + xlines) * 2 + 0) * 3 + 2] = 0;
      data[((i + xlines) * 2 + 1) * 3 + 0] = 0.5f;
      data[((i + xlines) * 2 + 1) * 3 + 1] = -0.5f + i / float(ylines - 1);
      data[((i + xlines) * 2 + 1) * 3 + 2] = 0;
    }

    return data;
  }
};
} // namespace primitives
} // namespace objects
} // namespace GL
} // namespace magnet
