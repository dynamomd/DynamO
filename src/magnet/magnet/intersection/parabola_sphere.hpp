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
#include <magnet/math/cubic.hpp>
#include <magnet/math/bisect.hpp>
#include <magnet/math/vector.hpp>

namespace magnet {
  namespace intersection {
    namespace detail {
      struct QuarticFunc
      {
      public:
	inline double operator()(double t) const
	{
	  return (((coeffs[0] * t + coeffs[1]) * t + coeffs[2]) * t + coeffs[3]) * t + coeffs[4];
	}
	
	double coeffs[5];
      };
    }

    /*! \brief A parabolic(ray)-sphere intersection test with backface culling.
      
      \param T The origin of the ray relative to the sphere center.
      \param D The direction/velocity of the ray.
      \param G The acceleration of the ray.
      \param r The radius of the sphere.
      \return The time until the intersection, or HUGE_VAL if no intersection.
    */
    inline double parabola_sphere_bfc(const math::Vector& T,
				      const math::Vector& D,
				      const math::Vector& G,
				      const double& r)
    {
      //This is our lengthscale to which we bisect the roots
      const double rootthreshold = 1e-16 * r;
      
      double DdotT = (D | T);

      detail::QuarticFunc f;
      f.coeffs[0] = 0.25 * G.nrm2();
      f.coeffs[1] = G | D;
      f.coeffs[2] = D.nrm2() + (G | T);
      f.coeffs[3] = 2 * DdotT;
      f.coeffs[4] = T.nrm2()- r * r;
      
      //We calculate the roots of the cubic differential of F
      //\f$F=A t^4 + B t^3 + C t^2 + D t + E == 0\f$ taking the differential gives
      //\f$F=4 A t^3 + 3 B t^2 + 2C t + D == 0\f$ and normalizing the cubic term gives
      //\f$F=t^3 + \frac{3 B}{4 A} t^2 + \frac{2C}{4 A} t + \frac{D}{4 A} == 0\f$
      double roots[3];
      size_t rootCount = magnet::math::cubicSolve(f.coeffs[1] * 3 / (4 * f.coeffs[0]),
						  f.coeffs[2] * 2 / (4 * f.coeffs[0]),
						  f.coeffs[3] * 1 / (4 * f.coeffs[0]),
						  roots[0], roots[1], roots[2]);
      
      //Sort the roots in ascending order
      std::sort(roots, roots + rootCount);
      if ((roots[0] > 0) && (f(roots[0]) < 0))
	{
	  if (f(0) <= 0) 
	    return 0;
	  else
	    return magnet::math::bisect(f, 0, roots[0], rootthreshold);
	}

      if ((rootCount == 3) && (roots[2] > 0) && (f(roots[2]) < 0))
	{
	  double tmin = std::max(0.0, roots[1]);
	  if (f(tmin) <= 0)
	    return tmin;
	  else
	    return magnet::math::bisect(f, tmin, roots[2], rootthreshold);
	}
      return HUGE_VAL;
    }

    /*! \brief A parabolic(ray) inverse-sphere intersection test with
        backface culling.
      
      \param T The origin of the ray relative to the inverse-sphere center.
      \param D The direction/velocity of the ray.
      \param G The acceleration of the ray.
      \param r The radius of the sphere.
      \return The time until the intersection, or HUGE_VAL if no intersection.
    */
    inline double parabola_invsphere_bfc(const math::Vector& T,
					 const math::Vector& D,
					 const math::Vector& G,
					 const double& r)
    {
      //This is our lengthscale to which we bisect the roots
      const double rootthreshold = 1e-16 * r;
      
      double DdotT = (D | T);

      detail::QuarticFunc f;
      f.coeffs[0] = -0.25 * G.nrm2();
      f.coeffs[1] = -G | D;
      f.coeffs[2] = -D.nrm2() + (G | T);
      f.coeffs[3] = -2 * DdotT;
      f.coeffs[4] = r * r - T.nrm2();
      
      //We calculate the roots of the cubic differential of F
      //\f$F=A t^4 + B t^3 + C t^2 + D t + E == 0\f$ taking the differential gives
      //\f$F=4 A t^3 + 3 B t^2 + 2C t + D == 0\f$ and normalizing the cubic term gives
      //\f$F=t^3 + \frac{3 B}{4 A} t^2 + \frac{2C}{4 A} t + \frac{D}{4 A} == 0\f$
      double roots[3];
      size_t rootCount = magnet::math::cubicSolve(f.coeffs[1] * 3 / (4 * f.coeffs[0]),
						  f.coeffs[2] * 2 / (4 * f.coeffs[0]),
						  f.coeffs[3] * 1 / (4 * f.coeffs[0]),
						  roots[0], roots[1], roots[2]);
      
      //Sort the roots in ascending order
      std::sort(roots, roots + rootCount);

      if ((rootCount == 3) && (roots[1] > 0) && (f(roots[1]) < 0))
	{
	  double tmin = std::max(0.0, roots[0]);
	  if (f(tmin) <= 0)
	    return tmin;
	  else
	    return magnet::math::bisect(f, tmin, roots[1], rootthreshold);
	}
      
      double tlast = roots[rootCount - 1];
      
      if (f(tlast) <= 0) return std::max(0.0, tlast);

      if ((tlast < 0) && (f(0) <= 0)) return 0;
      
      double t0 = std::max(0.0, tlast);
      
      double te = std::pow(4 * r * r / G.nrm2(), 0.25);
      
      while (f(te) >= 0)
	{
	  t0 += te;
	  te *= 2;
	}
      
      return magnet::math::bisect(f, t0, te, rootthreshold);
    }
  }
}
