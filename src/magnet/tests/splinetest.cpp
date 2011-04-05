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

#include <magnet/math/spline.hpp>
#include <iostream>
#include <fstream>

int main(int argc, char *argv[])
{
  magnet::math::Spline spline;

  spline.addPoint(0, 0);
  spline.addPoint(0.75, 0.5);
  spline.addPoint(1, 1);  

  {
    std::ofstream of("orig.dat");
    for (magnet::math::Spline::const_iterator iPtr = spline.begin();
	 iPtr != spline.end(); ++iPtr)
      of << iPtr->first << " " << iPtr->second << "\n";
  }
  
  {
    std::ofstream of("spline.natural.dat");
    for (double x(-0.05); x <= 1.06; x += 0.005)
      of << x << " " << spline(x) << "\n";
  }

  {
    std::ofstream of("spline.fixedy1.dat");
    spline.setBoundaryConditions(0, 0, magnet::math::Spline::FIXED_1ST_DERIV_BC);
    for (double x(-0.05); x <= 1.06; x += 0.005)
      of << x << " " << spline(x) << "\n";
  }

}
