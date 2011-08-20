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
#include <dynamo/dynamics/overlapFunc/CubePlane.hpp>
#include <boost/foreach.hpp>

bool 
dynamo::OverlapFunctions::CubePlane(const Vector& CubeOrigin, 
				    const Vector& CubeDimensions,
				    const Vector& PlaneOrigin, 
				    const Vector& PlaneNormal,
				    const double tol)
{
  Vector  relpos(CubeOrigin - PlaneOrigin);
  
  size_t counter[3] = {0, 0, 0};
  
  while (counter[NDIM-1] < 2)
    {
      Vector  pointpos(relpos);
      
      for (size_t iDim(0); iDim < NDIM; ++iDim)
	pointpos[iDim] += counter[iDim] * CubeDimensions[iDim];

      if ((pointpos | PlaneNormal) < tol) return true;
      
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
