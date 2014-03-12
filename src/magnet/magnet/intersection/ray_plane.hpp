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
#include <magnet/intersection/polynomial.hpp>

namespace magnet {
  namespace intersection {
    /*! \brief A stable particle and thick plane intersection test.
    
      \param R The origin of the particle relative to a location on the plane.
      \param V The direction/velocity of the particle.
      \param N The normal of the plane.
      \param d The interaction distance to the plane (the plane's thickness).
      \return The time until the intersection, or HUGE_VAL if no intersection.
    */
    inline double ray_plane(const math::Vector& R, const math::Vector& V, math::Vector N, const double d)
    {
      double r = R | N;
      if (r < 0) { r = -r; N = -N; }
      detail::PolynomialFunction<1> f(r - d, V | N);
      return detail::nextEvent(f);
    }
  }
}
