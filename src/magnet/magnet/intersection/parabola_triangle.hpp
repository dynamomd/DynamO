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
#include <magnet/intersection/parabola_plane.hpp>
#include <magnet/overlap/point_triangle.hpp>

namespace magnet {
  namespace intersection {
    /*! \brief A parabola-triangle intersection test.
      
      Here we make the assumption that all positions are relative to
      the first vertex. Thus we only need to pass the two edge vectors
      of the triangle.
      
      \param T The origin of the ray relative to the first vertex.
      \param D The direction/velocity of the ray.
      \param A The acceleration of the ray.
      \param E1 The first edge vector of the triangle (V1-V0).
      \param E2 The second edge vector of the triangle (V2-V0).
      \return The time until the intersection, or HUGE_VAL if no intersection.
    */
    inline double parabola_triangle_bfc(math::Vector T, 
					const math::Vector& D,
					const math::Vector& A,
					const math::Vector& E1, 
					const math::Vector& E2,
					const double d)
    {
      M_throw() << "Needs reimplementing now that parabola plane takes a diameter";
      math::Vector N = E1 ^ (E2 - E1);
      double t = magnet::intersection::parabola_plane_bfc(T, D, A, N, d);

      if (t == HUGE_VAL) return HUGE_VAL;

      T += t * D + 0.5 * t * t * A;
      
      if (magnet::overlap::point_triangle(T, E1, E2))
	return t;
      else
	return HUGE_VAL;
    }
  }
}
