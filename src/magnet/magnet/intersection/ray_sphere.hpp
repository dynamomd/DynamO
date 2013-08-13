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
    /*! \brief A ray-sphere intersection test.
       
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

      return std::max(0.0, - c / (TD - std::sqrt(arg)));
    }

    /*! \brief A ray-inverse_sphere intersection test.
      
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

    /*! \brief A ray-sphere intersection test where the sphere
      diameter is growing linearly with time.
      
      This uses only a slightly specialised form of the general stable
      EDMD algorithm, so the inverse sphere test is obtained by simply
      flipping the sign of the quadratic coefficients (through the
      inverse template parameter).

      \tparam inverse If true, this is the ray_inverse_growing_sphere test.
      \param R The origin of the ray relative to the sphere center.
      \param V The direction/velocity of the ray.
      \param r The radius of the sphere.
      \param inv_gamma The expansion rate of the sphere.
      \param t_curr The time passed since the sphere had a radius of r.
      \return The time until the intersection, or HUGE_VAL if no intersection.
    */
    template<bool inverse = false>
    inline double ray_growing_sphere(const math::Vector& R, const math::Vector& V, const double& r, const double inv_gamma, const double t_curr)
    {
      if (inv_gamma == 0)
	{
	  if (inverse) 
	    return ray_inv_sphere(R, V, r);
	  else
	    return ray_sphere(R, V, r);
	}

      const double fsign = (1 - 2 * inverse);
      const double a = fsign * ((V | V) - inv_gamma * inv_gamma * r * r);
      const double b = fsign * 2 * ((R | V) - r * r * inv_gamma * (1 + t_curr * inv_gamma));
      const double c = fsign * ((R | R) - (1 + t_curr * inv_gamma) * (1 + t_curr * inv_gamma) * r * r);
      const double arg = b * b - 4 * a * c;
      //Calculate the turning point.
      const double delta_t_min = - b / (2 * a);
      
      if (a == 0) return -c / b;
      
      //If a is negative, we will always have an interaction (the
      //overlap function will always go negative eventually) and it
      //will occur after the the turning point and after the latest
      //root if it exists. If a is positive, the interactions only
      //happen in the window between the first root and the turning
      //point.
      if (a < 0)
	{
	  //If there are no roots, then there are no intersections,
	  //but we will always interact once we pass the minimum.
	  if (arg < 0) return std::max(0.0, delta_t_min);

	  //We have roots, find the latest one
	  const double q = - b - std::copysign(std::sqrt(arg), b);
	  const double max_root = std::max(q / (2 * a), (2 * c) / q);
	  return std::max(max_root, 0.0);
	}
      else //a>0
	{
	  //If there are no roots, then there are no intersections,
	  if (arg < 0) return HUGE_VAL;
	  
	  const double q = - b - std::copysign(std::sqrt(arg), b);
	  const double min_root = std::min(q / (2 * a), (2 * c) / q);

	  if (delta_t_min < 0) return HUGE_VAL;
	  return std::max(min_root, 0.0);
	}
    }
  }
}
