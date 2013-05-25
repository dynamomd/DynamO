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
#include <dynamo/inputplugins/cells/cell.hpp>
#include <random>

namespace dynamo {
  struct CURandom: public UCell
  {
    CURandom(size_t nN, Vector  ndimensions, 
	     UCell* nextCell):
      UCell(nextCell),
      N(nN),
      dimensions(ndimensions),
      _rng(std::random_device()())
    {}

    size_t N;
    Vector dimensions;
    std::mt19937 _rng;

    virtual std::vector<Vector> placeObjects(const Vector & centre)
    {
      std::vector<Vector> retval;
      
      std::uniform_real_distribution<> uniform_dist;
      for (size_t i(0); i < N; ++i)
	{
	  Vector  position;
	  for (size_t iDim = 0; iDim < NDIM; iDim++)
	    position[iDim] = centre[iDim] - (uniform_dist(_rng) - 0.5) * dimensions[iDim];
	
	  //Get the next unit cells positions and push them to your list
	  const std::vector<Vector>& newsites = uc->placeObjects(position);
	  retval.insert(retval.end(), newsites.begin(), newsites.end());
	}

      return retval;    
    }
  };
}
