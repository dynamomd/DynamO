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
#include <magnet/collision/ray_triangle.hpp>

namespace magnet {
  namespace intersection {
    //! \brief A ray-quadrangle intersection test.
    //!
    //! The method used here is described in "An Efficient
    //! Ray-Quadrilateral Intersection Test", by Ares Lagae and Philip
    //! Dutré.
    //!
    //! \tparam BACKFACE_CULLING Ignores ray plane intersections
    //! where the ray enters the back of the plane.
    //! \param T The origin of the ray relative to a point on the plane.
    //! \param D The direction/velocity of the ray.
    //! \param N The normal of the plane.
    //! \param E1 An edge of the quadrangle (V1 - V0);
    //! \param E2 Another edge of the quadrangle (V2 - V0);
    //! \return The time until the intersection, or HUGE_VAL if no intersection.
    template<bool BACKFACE_CULLING>
    inline double ray_quadrilateral(const math::Vector& T,
				    const math::Vector& D,
				    const math::Vector& E1,
				    const math::Vector& E2)
    {
      //Do a test against the triangle of this quad, without the
      //diagonal test!
      double t1 = ray_triangle<BACKFACE_CULLING, false>(T, D, E1, E2);
      
      //Now test against the other triangle of the quad, also without
      //the diagonal test!
      double t2 = ray_triangle<BACKFACE_CULLING, false>(T - E1 - E2, D, -E1, -E2);
      
      //Just average the values, any HUGE_VAL's will propagate through
      //this operation.
      return (t1 + t2) * 0.5;
    }
  }
}
