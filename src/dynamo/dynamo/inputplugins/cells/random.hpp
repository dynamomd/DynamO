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
#include <dynamo/inputplugins/cells/cell.hpp>
#include <dynamo/base/is_simdata.hpp>

namespace dynamo {
  struct CURandom: public UCell
  {
    CURandom(size_t nN, Vector  ndimensions, 
	     boost::uniform_01<dynamo::baseRNG, double>& RNG,
	     UCell* nextCell):
      UCell(nextCell),
      N(nN),
      dimensions(ndimensions),
      uniform_sampler(RNG)
    {}

    size_t N;
    Vector  dimensions;
    boost::uniform_01<dynamo::baseRNG, double>& uniform_sampler;

    virtual std::vector<Vector  > placeObjects(const Vector & centre)
    {
      std::vector<Vector  > retval;

      for (size_t i(0); i < N; ++i)
	{
	  Vector  position;
	  for (size_t iDim = 0; iDim < NDIM; iDim++)
	    position[iDim] = centre[iDim] - (uniform_sampler() - 0.5) * dimensions[iDim];
	
	  //Get the next unit cells positions and push them to your list
	  const std::vector<Vector>& newsites = uc->placeObjects(position);
	  retval.insert(retval.end(), newsites.begin(), newsites.end());
	}

      return retval;    
    }
  };
}
