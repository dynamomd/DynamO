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
/*! \brief A point-triangle overlap test.

  This function assumes the point passed lies somewhere in the
  plane of the triangle and is defined relative to the first
  vertex of the triangle (V0).

  See http://softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm

  \param P The point's position, relative to V0.
  \param E1 The first edge vector of the triangle (V1-V0).
  \param E2 The second edge vector of the triangle (V2-V0).
  \return Whether the point is inside the triangle.
*/
inline bool point_triangle(const math::Vector &P, const math::Vector &E1,
                           const math::Vector &E2) {
  double uu = E1 | E1;
  double uv = E1 | E2;
  double vv = E2 | E2;
  double wu = P | E1;
  double wv = P | E2;
  double denom = uv * uv - uu * vv;

  double s = (uv * wv - vv * wu) / denom;
  if (s < 0.0 || s > 1.0)
    return false;

  double t = (uv * wu - uu * wv) / denom;
  if (t < 0.0 || (s + t) > 1.0)
    return false;

  return true;
}
} // namespace overlap
} // namespace magnet
