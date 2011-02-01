/*  DYNAMO:- Event driven molecular dynamics simulator 
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
    inline void MarcustoRGB(cl_uchar4& color, float h) 
    {
      h = clamp(h, 0.0f, 1.0f);

      float R = clamp(2.0f * h - 0.84f, 0.0f, 1.0f);
      float B = clamp(std::fabs(2.0f * h - 0.5f), 0.0f, 1.0f);

      float G;
      if (h < 0.3)
	G = 4.0 * h;
      else if (h < 0.92)
	G = 1.84 - (2.0 * h);
      else 
	G = (h / 0.08) - 11.5;

      G = clamp(G, 0.0f, 1.0f);

      color.s[0] = R * 255;
      color.s[1] = G * 255;
      color.s[2] = B * 255;
      color.s[3] = 255;
    }
  }
}
