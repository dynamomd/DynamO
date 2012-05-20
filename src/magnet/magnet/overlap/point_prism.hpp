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
#include <magnet/overlap/point_triangle.hpp>

namespace magnet {
  namespace overlap {
    //! \brief A point-prism (extruded triangle) overlap test.
    //!
    //! \param P The point's position, relative to V0.
    //! \param E1 The first edge vector of the prism face triangle (V1-V0).
    //! \param E2 The second edge vector of the prism face triangle (V2-V0).
    //! \param N Normal of the triangle.
    //! \param d The depth of the prism.
    //! \return Whether the point is inside the triangle.
    inline bool point_prism(const math::Vector& P, 
			    const math::Vector& E1, 
			    const math::Vector& E2,
			    const math::Vector& N,
			    const double d)
    {

      //Check the point is in the allowed depth range of the prism
      double dr = N | P;
      if ((dr > 0) || (dr < -d)) return false;

      return point_triangle(P, E1, E2);
    }
  }
}
