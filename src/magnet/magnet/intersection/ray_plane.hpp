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
    /*! \brief A ray-plane intersection test.
    
     \param T The origin of the ray relative to a point on the plane.
     \param D The direction/velocity of the ray.
     \param N The normal of the plane.
     \param d The interaction distance to the plane.
     \return The time until the intersection, or HUGE_VAL if no intersection.
    */
    inline double ray_plane(const math::Vector& T,
			    const math::Vector& D,
			    const math::Vector& N,
			    const double d)
    {
      double v = (D | N);
      double r = (T | N);
      if (//Particles must move to intersect
	  (v == 0)
	  //Particles must be moving towards each other to intersect
	  || ((v < 0) && (r < 0))
	  || ((v > 0) && (r > 0)))
	return HUGE_VAL;
      
      return  std::max(- (r - d) / v, 0.0);
    }
  }
}
