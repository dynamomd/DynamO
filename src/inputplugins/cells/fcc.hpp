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
#include "cell.hpp"

struct CUFCC: public CUCell
{
  CUFCC(CVector<long> ncells, Vector  ndimensions, CUCell* nextCell):
    CUCell(nextCell),
    cells(ncells),
    dimensions(ndimensions)
  {
    if (NDIM != 3) D_throw() << "FCC in other than 3 dims not allowed";
  }

  CVector<long> cells;
  Vector  dimensions;

  virtual std::vector<Vector  > placeObjects(const Vector & centre)
  {
    std::vector<Vector  > retval;

    Vector  cellWidth;
    
    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      cellWidth[iDim] = dimensions[iDim] / cells[iDim];
    
    Iflt rcoord[4][3];
    
    // sublattice a ---
    rcoord[0][0] = 0.0;
    rcoord[0][1] = 0.0;
    rcoord[0][2] = 0.0;
    
    // sublattice b ---
    rcoord[1][0] = cellWidth[0] * 0.5;
    rcoord[1][1] = cellWidth[1] * 0.5;
    rcoord[1][2] = 0.0;
    
    // sublattice c ---
    rcoord[2][0] = 0.0;
    rcoord[2][1] = cellWidth[1] * 0.5;
    rcoord[2][2] = cellWidth[2] * 0.5;
    
    // sublattice d ---
    rcoord[3][0] = cellWidth[0] * 0.5;
    rcoord[3][1] = 0.0;
    rcoord[3][2] = cellWidth[2] * 0.5;
    
    Vector  position;
    CVector<int> iterVec;
    
    for (iterVec[2] = 0; iterVec[2] < cells[2]; iterVec[2]++)
      for (iterVec[1] = 0; iterVec[1] < cells[1]; iterVec[1]++)
	for (iterVec[0] = 0; iterVec[0] < cells[0]; iterVec[0]++)
	  for (int iRef = 0; iRef < 4; iRef++)
	    {
	      for (size_t iDim = 0; iDim < NDIM; iDim++)
		position[iDim] =
		  rcoord[iRef][iDim] + cellWidth[iDim] * iterVec[iDim] - 0.5 * dimensions[iDim] + centre[iDim];
	      
	      //Get the next unit cells positions and push them to your list
	      BOOST_FOREACH(const Vector & vec, uc->placeObjects(position))
		retval.push_back(vec);
	    }

    return retval;    
  }
};
