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

namespace magnet {
  namespace intersection {
    /*! \brief A stable particle and thick plane intersection test.
    
      \param T The origin of the particle relative to a location on the plane.
      \param D The direction/velocity of the particle.
      \param N The normal of the plane.
      \param d The interaction distance to the plane (the plane's thickness).
      \return The time until the intersection, or HUGE_VAL if no intersection.
    */
    inline double ray_plane(const math::Vector& T, const math::Vector& D, const math::Vector& N, const double d)
    {
      double r = (T | N);
      double v = (D | N);

      if (r < 0) { r = -r; v = -v; }
      
      //Objects must move towards each other to intersect
      if (v >= 0) return HUGE_VAL;

      return std::max(-(r - d) / v, 0.0);
    }
  }
}
