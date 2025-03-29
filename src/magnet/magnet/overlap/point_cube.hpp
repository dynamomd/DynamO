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
/*! \brief Discovers if a cube and a point intersect.

  \param CubeOrigin The location of the cube relative to the point.
  \param CubeDimensions The size of the cube sides.
  \return If the point is inside/on the cube.
 */
inline bool point_cube(const magnet::math::Vector &CubeOrigin,
                       const magnet::math::Vector &CubeDimensions,
                       const double tol = 0) {
  for (size_t iDim(0); iDim < NDIM; ++iDim)
    if (fabs(CubeOrigin[iDim]) > (CubeDimensions[iDim] / 2))
      return false;

  return true;
}
} // namespace overlap
} // namespace magnet
