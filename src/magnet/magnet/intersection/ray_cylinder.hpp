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
    /*! \brief A ray_cylinder intersection test.
      
      \tparam inverse If true, the time the ray escapes the cylinder is returned.
      \param R The origin of the ray relative to a point on the cylinder axis.
      \param V The direction/velocity of the ray.
      \param A A normalized vector parallel to the axis of the cylinder.
      \param sig Radius of the cylinder.
            
      \return The time until the intersection, or HUGE_VAL if no intersection.
    */
    template<bool inverse = false>
    inline double ray_cylinder(math::Vector R, math::Vector V, const math::Vector& A, const double sig)
    {
      //Project off the axial component of the position and velocity
      R -= math::Vector((R | A) * A);
      V -= math::Vector((V | A) * A);
      return ray_sphere<inverse>(R, V, sig);
    }
  }
}
