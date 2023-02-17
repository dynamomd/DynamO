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
  struct CUHCP : public UCell
  {
    CUHCP(std::array<long, 3> ncells, Vector ndimensions, UCell *nextCell)
        : UCell(nextCell),
          cells(ncells)
    { 
      Vector _lattice_size{1.0, std::sqrt(3.0), 2.0 * std::sqrt(6.0) / 3.0};

      _cellDim = _lattice_size / std::max(_lattice_size[0], std::max(_lattice_size[1], _lattice_size[2]));
      _cellDim = elementwiseMultiply(_cellDim, ndimensions);
    }

    std::array<long, 3> cells;

    virtual std::vector<Vector> placeObjects(const Vector &centre)
    {
      Vector _lattice_size{1.0, std::sqrt(3.0), 2.0 * std::sqrt(6.0) / 3.0};
      double lattice_positions[4][3] = {
        {0,   0,                    0},
        {0.5, std::sqrt(3.0)/2,     0.0},
        {0.5, std::sqrt(3.0)/6,     std::sqrt(6.0) / 3.0},
        {0.0, 2.0 * std::sqrt(3.0) / 3, std::sqrt(6.0) / 3.0}
      };

      std::vector<Vector> retval;

      Vector cellWidth;
      for (size_t iDim = 0; iDim < NDIM; ++iDim)
        cellWidth[iDim] = _cellDim[iDim] / cells[iDim];

      Vector position;
      std::array<long, 3> iterVec = {{0, 0, 0}};

      while (iterVec[NDIM - 1] != cells[NDIM - 1])
      {
        //Create the zero position
        for (size_t iDim = 0; iDim < NDIM; iDim++)
          position[iDim] = cellWidth[iDim] * (static_cast<double>(iterVec[iDim]) /*+ 0.25 */) + centre[iDim];// - 0.5 * _cellDim[iDim];

        for (size_t i(0); i < 4; ++i) {
          Vector rel_pos{
            lattice_positions[i][0] / _lattice_size[0],
            lattice_positions[i][1] / _lattice_size[1],
            lattice_positions[i][2] / _lattice_size[2],
            };

          for (const Vector &vec : uc->placeObjects(position + elementwiseMultiply(rel_pos, cellWidth)))
            retval.push_back(vec);
        }

        // Now update the displacement vector
        iterVec[0]++;

        for (size_t iDim = 1; iDim < NDIM; iDim++)
        {
          // This increments the next dimension along when
          if (iterVec[iDim - 1] == cells[iDim - 1])
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
