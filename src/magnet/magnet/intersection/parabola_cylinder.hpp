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
#include <magnet/intersection/parabola_sphere.hpp>

namespace magnet {
namespace intersection {
/*! \brief A parabola-cylinder intersection test.

  This test ignores the back face of the cylinder, it is used to
  detect when a ray will enter a cylinder.

  \param T The origin of the ray relative to a point on the cylinder axis.
  \param D The direction/velocity of the ray.
  \param Aray The acceleration acting on the ray
  \param A A normalized vector parallel to the axis of the cylinder.
  \param r Radius of the cylinder.
  \return The time until the intersection, or HUGE_VAL if no intersection.
*/
inline double parabola_cylinder(math::Vector T, math::Vector D,
                                math::Vector Aray, const math::Vector &A,
                                const double r) {
  // Project off the axial component of the position, velocity and acceleration
  T -= math::Vector((T | A) * A);
  D -= math::Vector((D | A) * A);
  Aray -= math::Vector((Aray | A) * A);
  return parabola_sphere(T, D, Aray, r);
}
} // namespace intersection
} // namespace magnet
