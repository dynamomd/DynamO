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

#include "spline.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
  magnet::math::Spline spline;

  spline.addPoint(0,0);
  spline.addPoint(0.5,0.5);
  spline.addPoint(1,1);  
  spline.generate();

  std::cout << "\nOriginal points";
  for (std::vector<std::pair<double, double> >::const_iterator iPtr = spline.begin();
       iPtr != spline.end(); ++iPtr)
    std::cout << "\n" << iPtr->first << " " << iPtr->second;
  
  std::cout << "\nSplined points";
  for (double x(-0.05); x <= 1.06; x += 0.005)
    std::cout << "\n" << x << " " << spline(x);
}
