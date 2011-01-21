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
#include "cell.hpp"
#include "../../datatypes/vector.hpp"
#include <cmath>

struct CUMirror: public CUCell
{
  CUMirror(double F, CUCell* nextCell):
    CUCell(nextCell),
    fraction(F),
    count1(0),
    count2(0)
  {}

  ~CUMirror()
  {
    std::cout << "\nACTUAL CHIRALITY = " << (static_cast<double>(count1) / static_cast<double>(count1+count2));
  }

  double fraction;
  long count1, count2;

  virtual std::vector<Vector  > placeObjects(const Vector & centre)
  {
    //Must be placed at zero for the mirroring to work correctly
    std::vector<Vector  > retval(uc->placeObjects(Vector (0,0,0)));

    //Avoid dividing by zero, then distribute the images according to the fraction
    if (!(count1+count2) || (static_cast<double>(count1) / static_cast<double>(count1+count2) > fraction))
      ++count2;       
    else
      {
	++count1;
	//Mirror the unit cell now
	BOOST_FOREACH(Vector & vec, retval)
	  {
	    if (NDIM % 2)
	      //Odd dimensions, flip all for symmetry
		vec *= -1.0;
	    else
	      //Just exchange one dimension
	      vec[0] *= -1.0;	    
	  }
      }

    //Adjust the centering of the unit cell
    BOOST_FOREACH(Vector & vec, retval)
      vec += centre;

    return retval;    
  }
};
