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
#include <math.h>

namespace magnet {
  namespace intersection {
    /*! \brief A ray-sphere intersection test with backface culling.
       
      \param T The origin of the ray relative to the sphere center.
      \param D The direction/velocity of the ray.
      \param r The radius of the sphere.
      \return The time until the intersection, or HUGE_VAL if no intersection.
    */
    inline double ray_sphere(const math::Vector& T, const math::Vector& D, const double& r)
    {
      double TD = (T | D);

      if (TD >= 0) return HUGE_VAL;
      
      double c = T.nrm2() - r * r;
      double arg = TD * TD - D.nrm2() * c;
      
      if (arg < 0) return HUGE_VAL;

      return  std::max(0.0, - c / (TD - std::sqrt(arg)));
    }

    /*! \brief A ray-inverse_sphere intersection test with backface culling.
      
      An inverse sphere means an "enclosing" sphere.
      
      \param T The origin of the ray relative to the inverse sphere
      center.
      \param D The direction/velocity of the ray.
      \param d The diameter of the inverse sphere.

      \return The time until the intersection, or HUGE_VAL if no intersection.
    */
    inline double ray_inv_sphere(const math::Vector& T, const math::Vector& D, const double& r)
    {
      double D2 = D.nrm2();      
      double TD = T | D;
      double c = r * r - T.nrm2();
      double arg = TD * TD + D2 * c;

      if (D2 == 0) return HUGE_VAL;

      if (arg >= 0)
	{
	  double q = TD + copysign(std::sqrt(arg), TD);
	  return std::max(0.0, std::max(- q / D2, c / q));
	}
	
      //The ray never passes through the sphere, return the time when
      //it is closest to the sphere surface
      return std::max(0.0, - TD / D2);
    }
  }
}
