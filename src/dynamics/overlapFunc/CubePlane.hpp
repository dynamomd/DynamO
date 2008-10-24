/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "../../datatypes/vector.hpp"
#include <boost/foreach.hpp>

namespace DYNAMO
{
  namespace OverlapFunctions
  {
    /*! \brief Discovers if a cube and a plane intersect by testing
     * which side of the plane the points of the cube lie on.
     *
     * This is used in the collision CLSentinel to install itself in cells
     */
    bool CubePlane(const CVector<>& CubeOrigin, const CVector<>& CubeDimensions,
		   const CVector<>& PlaneOrigin, const CVector<>& PlaneNormal)
    {
      CVector<> relpos(CubeOrigin - PlaneOrigin);
      
      //Get which side the cube origin is on 
      bool OriginSign(std::signbit(relpos % PlaneNormal));

      CVector<size_t> counter(0);
      
      while (counter[NDIM-1] < 2)
	{
	  CVector<> pointpos(relpos);
	  
	  for (size_t iDim(0); iDim < NDIM; ++NDIM)
	    if (counter[iDim]) pointpos[iDim] += CubeDimensions[iDim];

	  if (std::signbit(pointpos % PlaneNormal) != OriginSign)
	    return true;
	  
	  ++counter[0];

	  for (size_t iDim(0); iDim < NDIM-1; ++iDim)
	    if (counter[iDim] > 1)
	      {
		counter[iDim] = 0;
		++counter[iDim+1];
	      }
	}

      return false;
    }
  };
};
