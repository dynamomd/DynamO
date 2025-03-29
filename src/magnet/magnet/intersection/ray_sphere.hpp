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
#include <magnet/intersection/polynomial.hpp>
#include <magnet/math/vector.hpp>

namespace magnet {
namespace intersection {
/*! \brief A ray-sphere intersection test.

  \tparam inverse If true, this returns the time the ray escapes the sphere
  (rather than enters).
  \param R The origin of the ray relative to the sphere center.
  \param V The direction/velocity of the ray.
  \param sig The radius of the sphere.
  \return The time until the intersection, or HUGE_VAL if no intersection.
*/
template <bool inverse = false>
inline double ray_sphere(const math::Vector &R, const math::Vector &V,
                         const double &sig) {
  detail::PolynomialFunction<2> f(R.nrm2() - sig * sig, 2 * (R | V),
                                  2 * V.nrm2());
  if (inverse)
    f.flipSign();
  return detail::nextEvent(f);
}

/*! \brief A ray-sphere intersection test where the sphere
  diameter is growing linearly with time.

  This uses only a slightly specialised form of the general stable
  EDMD algorithm, so the inverse sphere test is obtained by simply
  flipping the sign of the quadratic coefficients (through the
  inverse template parameter).

  \tparam inverse If true, this is the ray_inverse_growing_sphere test.
  \param R The origin of the ray relative to the sphere center.
  \param V The direction/velocity of the ray.
  \param sig The radius of the sphere.
  \param inv_gamma The expansion rate of the sphere.
  \param t_curr The time passed since the sphere had a radius of r.
  \return The time until the intersection, or HUGE_VAL if no intersection.
*/
template <bool inverse = false>
inline double ray_growing_sphere(const math::Vector &R, const math::Vector &V,
                                 const double &sig, const double inv_gamma,
                                 const double t_curr) {
  const double currentDiam = sig * (1 + inv_gamma * t_curr);
  detail::PolynomialFunction<2> f(
      R.nrm2() - currentDiam * currentDiam,
      2 * (R | V) - 2 * inv_gamma * sig * currentDiam,
      2 * (V.nrm2() - sig * sig * inv_gamma * inv_gamma));
  /*We cannot determine the sign of f2 at compile time*/
  if (inverse)
    f.flipSign();
  return detail::nextEvent(f);
}
} // namespace intersection
} // namespace magnet
