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
    //! \brief A ray-plane intersection test.
    //!
    //! \tparam BACKFACE_CULLING Ignores ray plane intersections
    //! where the ray enters the back of the plane.
    //! \param T The origin of the ray relative to a point on the plane.
    //! \param D The direction/velocity of the ray.
    //! \param N The normal of the plane.
    //! \return The time until the intersection, or HUGE_VAL if no intersection.
    template<bool BACKFACE_CULLING>
    inline double ray_plane(const math::Vector& T,
			    const math::Vector& D,
			    const math::Vector& N)
    {
      double ND = (D | N);
      
      if (BACKFACE_CULLING)
	{ if (ND >= 0) return HUGE_VAL; }
      else
	{ if (ND == 0) return HUGE_VAL; }

      return  -(T | N) / ND;
    }
  }
}
