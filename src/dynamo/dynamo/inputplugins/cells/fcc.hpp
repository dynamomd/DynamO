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

namespace dynamo
{
  struct CUFCC : public UCell
  {
    CUFCC(std::array<long, 3> ncells, Vector ndimensions, UCell *nextCell)
        : UCell(nextCell),
          cells(ncells)
    {
      _cellDim = ndimensions;
      if (NDIM != 3)
        M_throw() << "FCC in other than 3 dims not allowed";
    }

    std::array<long, 3> cells;

    virtual std::vector<Vector> placeObjects(const Vector &centre)
    {
      std::vector<Vector> retval;

      Vector cellWidth;

      for (size_t iDim = 0; iDim < NDIM; ++iDim)
        cellWidth[iDim] = _cellDim[iDim] / cells[iDim];

      double rcoord[4][3];

      // The unit cell is centered in the cell box
      //  sublattice a ---
      rcoord[0][0] = cellWidth[0] * 0.25;
      rcoord[0][1] = cellWidth[1] * 0.25;
      rcoord[0][2] = cellWidth[2] * 0.25;

      // sublattice b ---
      rcoord[1][0] = cellWidth[0] * 0.75;
      rcoord[1][1] = cellWidth[1] * 0.75;
      rcoord[1][2] = cellWidth[2] * 0.25;

      // sublattice c ---
      rcoord[2][0] = cellWidth[0] * 0.25;
      rcoord[2][1] = cellWidth[1] * 0.75;
      rcoord[2][2] = cellWidth[2] * 0.75;

      // sublattice d ---
      rcoord[3][0] = cellWidth[0] * 0.75;
      rcoord[3][1] = cellWidth[1] * 0.25;
      rcoord[3][2] = cellWidth[2] * 0.75;

      Vector position;
      std::array<int, 3> iterVec;

      for (iterVec[2] = 0; iterVec[2] < cells[2]; iterVec[2]++)
        for (iterVec[1] = 0; iterVec[1] < cells[1]; iterVec[1]++)
          for (iterVec[0] = 0; iterVec[0] < cells[0]; iterVec[0]++)
            for (int iRef = 0; iRef < 4; iRef++)
            {
              for (size_t iDim = 0; iDim < NDIM; iDim++)
                position[iDim] =
                    rcoord[iRef][iDim] + cellWidth[iDim] * iterVec[iDim] - 0.5 * _cellDim[iDim] + centre[iDim];

              // Get the next unit cells positions and push them to your list
              const std::vector<Vector> &newsites = uc->placeObjects(position);
              retval.insert(retval.end(), newsites.begin(), newsites.end());
            }

      return retval;
    }

  };
}
