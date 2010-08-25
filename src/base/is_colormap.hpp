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

#include <cmath>

namespace DYNAMO
{
  struct RGB {
    RGB(Iflt a, Iflt b, Iflt c): R(a),G(b),B(c) {}
    Iflt R;
    Iflt G;
    Iflt B;
  };

  template<class T>
  class ColorMap
  {
  public:
    ColorMap(T nstart, T nend):
      start(nstart),
      end(nend)
    {}
    
    RGB getColor(T val)
    {
      if ((end-start == 0.0) || (val < start))
	return RGB(0, 0, 1);
      else if (val > end)
	return RGB(1, 0, 0);
	      
      Iflt normval = static_cast<Iflt>(val - start) / static_cast<Iflt>(end - start);      

      Iflt R = 2.0 * normval - 0.84;
      Iflt G;
      Iflt B = std::fabs(2.0 * normval - 0.5);

      if (val < 0.3)
	G = 4.0 * normval;
      else if (normval < 0.92)
	G = 1.84 - (2.0 * normval);
      else 
	G = (normval / 0.08) - 11.5;

      boundVal(R);
      boundVal(G);
      boundVal(B);

      return RGB(R,G,B);
    }
    
  private:

    inline void boundVal(Iflt& val)
    {
      if (val > 1.0) 
	val = 1.0;
      else if (val < 0.0)
	val = 0.0;
    }

    T start;
    T end;
  };
}
