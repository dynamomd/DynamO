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
#include <magnet/intersection/ray_sphere.hpp>

namespace magnet {
  namespace intersection {
    /*! \brief A ray-inverse_cylinder intersection test.
    
     This test ignores the back face of the cylinder, it is used to
     detect when a ray inside a cylinder will escape.
     
     \param T The origin of the ray relative to a point on the cylinder axis.
     \param D The direction/velocity of the ray.
     \param A A normalized vector parallel to the axis of the cylinder.
     \param r Radius of the cylinder.

      \tparam always_intersect If true, this will ensure that glancing
      ray's never escape the enclosing cylinder by returning the time
      when the ray is nearest the cylinder if the ray does not intersect
      the cylinder.

     \return The time until the intersection, or HUGE_VAL if no intersection.
    */
    inline double ray_inv_cylinder_bfc(math::Vector T, 
				       math::Vector D, 
				       const math::Vector& A, 
				       const double r)
    {
      //Project off the axial component of the position and velocity
      T -= math::Vector((T | A) * A);
      D -= math::Vector((D | A) * A);

      return ray_inv_sphere_bfc(T, D, r);
    }


    //! \brief A ray-cylinder intersection test.
    //!
    //! This test ignores the back face of the cylinder, it is used to
    //! detect when a ray will enter a cylinder.
    //! 
    //! \param T The origin of the ray relative to a point on the cylinder axis.
    //! \param D The direction/velocity of the ray.
    //! \param A A normalized vector parallel to the axis of the cylinder.
    //! \param r Radius of the cylinder.
    //! \return The time until the intersection, or HUGE_VAL if no intersection.
    inline double ray_cylinder_bfc(math::Vector T, math::Vector D, const math::Vector& A, const double r)
    {
      //Project off the axial component of the position and velocity
      T -= math::Vector((T | A) * A);
      D -= math::Vector((D | A) * A);

      return ray_sphere_bfc(T, D, r);
    }
  }
}
