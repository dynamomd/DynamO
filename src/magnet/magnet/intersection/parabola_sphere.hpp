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
#include <magnet/math/cubic.hpp>
#include <magnet/math/bisect.hpp>
#include <magnet/math/vector.hpp>

namespace magnet {
  namespace intersection {
    namespace detail {
      struct QuarticFunc
      {
      public:
	inline double operator()(double t)
	{
	  return (((coeffs[0] * t + coeffs[1]) * t + coeffs[2]) * t + coeffs[3]) * t + coeffs[4];
	}
	
	double coeffs[5];
      };
    }

    //! A parabolic(ray)-sphere intersection test with backface culling.
    //!
    //! \param T The origin of the ray relative to the sphere center.
    //! \param D The direction/velocity of the ray.
    //! \param G The acceleration of the ray.
    //! \param r The radius of the sphere.
    //! \return The time until the intersection, or HUGE_VAL if no intersection.
    inline double parabola_sphere_bfc(const Vector& T,
				      const Vector& D,
				      const Vector& G,
				      const double& r)
    {
      //This is our lengthscale to which we bisect the roots
      const double rootthreshold = 1e-16 * r;
      
      double DdotT = (D | T);

      magnet::math::Bisect<detail::QuarticFunc> quartic;
      quartic.coeffs[0] = 0.25 * G.nrm2();
      quartic.coeffs[1] = G | D;
      quartic.coeffs[2] = D.nrm2() + (G | T);
      quartic.coeffs[3] = 2 * DdotT;
      quartic.coeffs[4] = T.nrm2()- r * r;
      
      //We calculate the roots of the cubic differential of F
      //\f$F=A t^4 + B t^3 + C t^2 + D t + E == 0\f$ taking the differential gives
      //\f$F=4 A t^3 + 3 B t^2 + 2C t + D == 0\f$ and normalizing the cubic term gives
      //\f$F=t^3 + \frac{3 B}{4 A} t^2 + \frac{2C}{4 A} t + \frac{D}{4 A} == 0\f$
      double roots[3];
      size_t rootCount = magnet::math::cubicSolve(quartic.coeffs[1] * 3 / (4 * quartic.coeffs[0]),
						  quartic.coeffs[2] * 2 / (4 * quartic.coeffs[0]),
						  quartic.coeffs[3] * 1 / (4 * quartic.coeffs[0]),
						  roots[0], roots[1], roots[2]);
      
      //Sort the roots in ascending order
      std::sort(roots, roots + rootCount);
      //The roots are (in order) a minimum, (and if we have 3 roots), a
      //local maximum, then another minimum
      
      //We need to define the dynamics of overlapping particles
      //
      //We follow the convention that, if particles are overlapped and
      //approaching, we collide them now.
      //
      //But, there are some cases where the particle is overlapped,
      //receeding but will not escape the overlap before turning around
      //and approaching again.
      //
      //Here we make another choice, we collide whenever the particle is
      //touching/overlapping and approaching. We may get some kind of
      //inelastic collapse, which is fair, but it may be avoided by
      //sleeping the particles or by other dynamics
      if (quartic(0) <= 0)
	{//We're overlapping! We need the overlapped dynamics!
	  if (DdotT < 0)
	    //We're overlapping and approaching! Cause an instantaneous collision
	    return 0;
      
	  //We're overlapping but receeding. If there's only one cubic
	  //root, then we're fine, we've passed the local minimum (we're
	  //receeding) so we'll escape the particle (never to return) soon
	  //enough, return no collision
	  if (rootCount == 1) return HUGE_VAL;

	  //We have a local minimum, maximum and minimum. We must be after
	  //the first local minimum as we're overlapping but receeding. If
	  //the maximum is in the past, and we're receeding, we must be
	  //past the second minimum and we will only get further away! so
	  //return no collision
	  if (roots[1] < 0) return HUGE_VAL;

	  //The local maximum is in the future (or now). Check if we will be
	  //outside the particle then.
	  if (quartic(roots[1]) > 0)
	    { //We will be outside the particle, so check if we re-enter it!
	      if (quartic(roots[2]) < 0)
		{//We do reenter! So return the time of this reentry as the next event
		  return std::max(0.0, quartic.bisectRoot(roots[1], roots[2], rootthreshold));
		}
	      //We wont reenter, so return no intersection
	      return HUGE_VAL;
	    }
      
	  //Our local maxima is in the future, and it's still an
	  //overlapped state!  Return an event at the local maximum, it's
	  //the furthest away we're going to get without approaching
	  //again. Either we'll have another event before we reach it or
	  //we'll have to handle it by sleeping the particle.
	  return roots[1]; //return the local maximum 
	}
      
      //We're not overlapping, so check for standard events  

      //Check the first minimum (we always have one), if it's in the
      //future check that we overlap then, and then bisect the event
      if ((roots[0] >= 0) && (quartic(roots[0]) <= 0))
	return std::max(0.0, quartic.bisectRoot(0, roots[0], rootthreshold));

      //Check the second minimum if we have one. When bisecting we use
      //either now or the local maximum, depending on which is later as
      //the first point.
      if ((rootCount > 1)
	  && (roots[2] > 0)
	  && (quartic(roots[2]) < 0)
	  && (quartic(std::max(0.0, roots[1])) >= 0))
	return std::max(0.0, quartic.bisectRoot(std::max(0.0, roots[1]), roots[2], rootthreshold));
  
      return HUGE_VAL;
    }
  }
}
