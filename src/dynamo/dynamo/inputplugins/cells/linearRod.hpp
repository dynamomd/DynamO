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
#include "cell.hpp"
#include "../../datatypes/vector.hpp"
#include <cmath>

struct CUlinearRod: public CUCell
{
  CUlinearRod(size_t pcl, double WL, CUCell* nextCell):
    CUCell(nextCell),
    pairchainlength(pcl),
    walklength(WL)
  {
    if (pcl == 0) M_throw() << "Cant have zero chain length";
  }

  size_t pairchainlength;  
  double walklength;
  
  virtual std::vector<Vector> placeObjects(const Vector & centre)
  {
    Vector tmp(0,0,0);

    std::vector<Vector> retval;

    for (size_t iStep = 0; iStep < pairchainlength; ++iStep)
      { 
	tmp[0] = (double(iStep) - (double(walklength) * 0.5)) * walklength;
	
	BOOST_FOREACH(const Vector & vec, uc->placeObjects(tmp + centre))
	  retval.push_back(vec);
      }

    return retval;
  }
};
