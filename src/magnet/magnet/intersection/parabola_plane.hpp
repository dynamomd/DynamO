/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
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

namespace magnet {
  namespace intersection {
    /*! \brief A parabola-plane intersection test which ignores
     * negative time intersections.
     *
     * \param T The origin of the ray relative to a point on the plane.
     * \param D The direction/velocity of the ray.
     * \param A The acceleration of the ray.
     * \param N The normal of the plane.
     * \return The time until the intersection, or HUGE_VAL if no intersection.
     */
    inline double parabola_plane_bfc(const Vector& T,
				     const Vector& D,
				     const Vector& A,
				     const Vector& N)
    {
      double adot = (N | A);
      if (adot == 0) return ray_plane<true>(T, D, N);

      double rdot = T | N;
      double vdot = D | N;
      double arg = vdot * vdot - 2 * rdot * adot;
      
      //Check if the particle curves away from the wall before it hits the wall
      if (arg < 0) return HUGE_VAL;

      //The particle will hit the wall

      double t = -(vdot + ((vdot<0) ? -1: 1) * std::sqrt(arg));
      double x1 = t / adot;
      double x2 = 2 * rdot / t;
      
      //If (adot > 0), the top of the arc of the particle passes
      //through the plate, so we want the earliest root. 
      //
      //If (adot < 0), the top of the arc of the particle is outside
      //the plate, then we want the earliest root

      return (adot > 0) ? std::min(x1, x2) : std::max(x1, x2);
    }
  }
}
