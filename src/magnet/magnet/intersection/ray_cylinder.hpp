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
#include <magnet/math/vector.hpp>

namespace magnet {
  namespace intersection {
    //! A ray-cylinder intersection test.
    //!
    //! \param T The origin of the ray relative to a point on the cylinder axis.
    //! \param D The direction/velocity of the ray.
    //! \param A A normalized vector parallel to the axis of the cylinder.
    //! \param r Radius of the cylinder.
    //! \return The time until the intersection, or HUGE_VAL if no intersection.
    inline double ray_cylinder(Vector T,
			       Vector D,
			       const Vector& A,
			       const double r)
    {
      T -= Vector((T | A) * A);
      D -= Vector((D | A) * A);

      double TD = T | D;
      double D2 = D.nrm2();
      double arg = TD * TD - D2 * (T.nrm2() - r * r);
      
      if (arg < 0) return HUGE_VAL;

      return  (std::sqrt(arg) - TD) / D2;
    }
  }
}
