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
#include "cell.hpp"
#include "../../datatypes/vector.hpp"
#include <cmath>

struct CUringSnake: public CUCell
{
  CUringSnake(size_t pcl, Iflt WL, CUCell* nextCell):
    CUCell(nextCell),
    pairchainlength(pcl),
    walklength(WL)
  {
    if (pcl == 0) M_throw() << "Cant have zero chain length";
  }

  size_t pairchainlength;  
  Iflt walklength;
  
  virtual std::vector<Vector  > placeObjects(const Vector & centre)
  {
    size_t L(size_t(std::sqrt(pairchainlength)));
    
    std::vector<Vector  > localsites;

    Vector  x(0,0,0);

    Iflt direction(walklength);
    
    for (size_t i(0); i < pairchainlength;)
      {
	if (i % L)
	  x[0] += direction;
	else
	  {
	    x[1] += walklength;
	    direction *= -1;
	  }

	localsites.push_back(x);
	
	++i;
      }

    direction *= -1;
    x[2] += walklength;

    for (size_t i(pairchainlength - 1); i != 0;)
      {
	localsites.push_back(x);
	
	if (i % L)
	  x[0] += direction;
	else
	  {
	    x[1] -= walklength;
	    direction *= -1;
	  }

	--i;
      }

    localsites.push_back(x);
        
    std::vector<Vector  > retval;
    BOOST_FOREACH(const Vector & vec, localsites)
      BOOST_FOREACH(const Vector & vec2, uc->placeObjects(vec))
        retval.push_back(vec2);

    return retval;    
  }
};
