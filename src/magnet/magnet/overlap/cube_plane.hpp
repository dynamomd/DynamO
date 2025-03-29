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
#include <magnet/math/vector.hpp>

namespace magnet {
namespace overlap {
/*! \brief Discovers if a cube and a plane intersect by testing
  which side of the plane the points of the cube lie on.

  This is used in the collision CLSentinel to install itself in
  cells
 */
inline bool cube_plane(const magnet::math::Vector &CubeOrigin,
                       const magnet::math::Vector &CubeDimensions,
                       const magnet::math::Vector &PlaneOrigin,
                       const magnet::math::Vector &PlaneNormal,
                       const double tol = 0) {
  magnet::math::Vector relpos(CubeOrigin - PlaneOrigin);

  size_t counter[3] = {0, 0, 0};

  while (counter[NDIM - 1] < 2) {
    magnet::math::Vector pointpos(relpos);

    for (size_t iDim(0); iDim < NDIM; ++iDim)
      pointpos[iDim] += counter[iDim] * CubeDimensions[iDim];

    if ((pointpos | PlaneNormal) < tol)
      return true;

    ++counter[0];

    for (size_t iDim(0); iDim < NDIM - 1; ++iDim)
      if (counter[iDim] > 1) {
        counter[iDim] = 0;
        ++counter[iDim + 1];
      }
  }

  return false;
}
} // namespace overlap
} // namespace magnet
