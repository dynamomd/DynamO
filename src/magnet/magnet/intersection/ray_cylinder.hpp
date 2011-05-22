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
#include <magnet/intersection/ray_sphere.hpp>

namespace magnet {
  namespace intersection {
    //! A ray-inverse_cylinder intersection test.
    //!
    //! This test ignores the back face of the cylinder, it is used to
    //! detect when a ray inside a cylinder will escape.
    //! 
    //! \param T The origin of the ray relative to a point on the cylinder axis.
    //! \param D The direction/velocity of the ray.
    //! \param A A normalized vector parallel to the axis of the cylinder.
    //! \param r Radius of the cylinder.
    //! \return The time until the intersection, or HUGE_VAL if no intersection.
    inline double ray_inv_cylinder_bfc(Vector T, 
				       Vector D, 
				       const Vector& A, 
				       const double r)
    {
      //Project off the axial component of the position and velocity
      T -= Vector((T | A) * A);
      D -= Vector((D | A) * A);

      return ray_inv_sphere(T, D, r);
    }


    //! A ray-inverse_cylinder intersection test.
    //!
    //! This test ignores the back face of the cylinder, it is used to
    //! detect when a ray will enter a cylinder.
    //! 
    //! \param T The origin of the ray relative to a point on the cylinder axis.
    //! \param D The direction/velocity of the ray.
    //! \param A A normalized vector parallel to the axis of the cylinder.
    //! \param r Radius of the cylinder.
    //! \return The time until the intersection, or HUGE_VAL if no intersection.
    inline double ray_cylinder_bfc(Vector T, Vector D, const Vector& A, const double r)
    {
      //Project off the axial component of the position and velocity
      T -= Vector((T | A) * A);
      D -= Vector((D | A) * A);

      return ray_sphere(T, D, r);
    }
  }
}
