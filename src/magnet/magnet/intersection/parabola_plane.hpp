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
#include <magnet/intersection/ray_plane.hpp>
#include <magnet/math/quadratic.hpp>

namespace magnet {
namespace intersection {
/*! \brief A parabola-plane intersection test which ignores
  negative time intersections.

  \param T The origin of the ray relative to a point on the plane.
  \param D The direction/velocity of the ray.
  \param A The acceleration of the ray.
  \param N The normal of the plane.
  \param d The thickness of the plane.
  \return The time until the intersection, or HUGE_VAL if no intersection.
 */
inline double parabola_plane(const math::Vector &R, const math::Vector &V,
                             const math::Vector &A, math::Vector N,
                             const double d) {
  double rdotn = N | R;
  if (rdotn < 0) {
    N = -N;
    rdotn = -rdotn;
  }
  detail::PolynomialFunction<2> f(rdotn - d, V | N, A | N);
  return detail::nextEvent(f);
}
} // namespace intersection
} // namespace magnet
