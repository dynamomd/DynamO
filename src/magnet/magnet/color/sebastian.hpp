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

#include <cmath>
#include <magnet/clamp.hpp>
#include <magnet/exception.hpp>

namespace magnet {
  namespace color {
    inline void SebastiantoRGB(GLfloat color[4], float h) 
    {
      h = clamp(h, 0.0f, 1.0f);

      color[0] = h;
      color[1] = 4 * h * (1 - h);
      color[2] = (1 - h);
      color[3] = 1;
    }
  }
}
