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

namespace dynamo {
  struct CUBCC: public UCell
  {
    CUBCC(std::tr1::array<long, 3> ncells, Vector  ndimensions, UCell* nextCell):
      UCell(nextCell),
      cells(ncells),
      dimensions(ndimensions)
    {}

    std::tr1::array<long, 3> cells;
    Vector  dimensions;

    virtual std::vector<Vector> placeObjects(const Vector & centre)
    {
      std::vector<Vector> retval;

      Vector  cellWidth;
      for (size_t iDim = 0; iDim < NDIM; ++iDim)
	cellWidth[iDim] = dimensions[iDim] / cells[iDim];
    
      Vector  position;
      std::tr1::array<long, 3> iterVec = {{0,0,0}};
    
      while (iterVec[NDIM - 1] != cells[NDIM-1])
	{      
	  for (size_t iDim = 0; iDim < NDIM; iDim++)
	    position[iDim] = cellWidth[iDim] * (static_cast<double>(iterVec[iDim]) + 0.25) - 0.5 * dimensions[iDim] 
	      + centre[iDim];
	
	  BOOST_FOREACH(const Vector & vec, uc->placeObjects(position))
	    retval.push_back(vec);
	
	  for (size_t iDim = 0; iDim < NDIM; iDim++)
	    position[iDim] += cellWidth[iDim]/2.0;
	
	  const std::vector<Vector>& newsites = uc->placeObjects(position);
	  retval.insert(retval.end(), newsites.begin(), newsites.end());
	
	  //Now update the displacement vector
	  iterVec[0]++;
	
	  for (size_t iDim = 1; iDim < NDIM; iDim++)
	    {
	      //This increments the next dimension along when
	      if (iterVec[iDim - 1] == cells[iDim -1])
		{
		  iterVec[iDim - 1] = 0;
		  iterVec[iDim] += 1;
		}
	    }
	}
    
      return retval;
    }
  };
}
