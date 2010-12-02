/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include <magnet/clamp.hpp>
#include <magnet/exception.hpp>

namespace magnet {
  namespace color {
    inline void HSVtoRGB(double& r, double& g, double& b, double h, double s = 1, double v = 1) 
    {
      double temp;
      h = modf(h, &temp);

      s = clamp(s, 0.0, 1.0);
      v = clamp(v, 0.0, 1.0);

      h = h * 6;

      unsigned int i = h;
      double f = h - i;
      double p = v * (1 - s);
      double q = v * (1 - s * f);
      double t = v * (1 - s * (1 - f));

      switch(i) {
      case 0:	
	r = v;
	g = t;
	b = p;
	break;

      case 1:	
	r = q;
	g = v;
	b = p;
	break;

      case 2:
	r = p;
	g = v;
	b = t;
	break;

      case 3:
	r = p;
	g = q;
	b = v;
	break;

      case 4:
	r = t;
	g = p;
	b = v;
	break;

      case 5:
	r = v;
	g = p;
	b = q;
	break;

      default:
	M_throw() << "Out of range H? Can not happen";
      }
    }
  }
}
