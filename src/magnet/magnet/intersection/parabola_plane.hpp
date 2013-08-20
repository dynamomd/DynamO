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
#include <magnet/intersection/ray_plane.hpp>
#include <magnet/math/quadratic.hpp>

namespace magnet {
  namespace intersection {
    /*! \brief A parabola-plane intersection test which ignores
      negative time intersections.
     
      \param T The origin of the ray relative to a point on the plane.
      \param D The direction/velocity of the ray.
      \param A The acceleration of the ray.
      \param N The normal of the plane.
      \param d The thickness of the plane.
      \return The time until the intersection, or HUGE_VAL if no intersection.
     */
    inline double parabola_plane(const math::Vector& T, const math::Vector& D, const math::Vector& A, math::Vector N, const double d)
    {
      //Check if this is actually a ray-plane test
      double adot = A | N;
      if (adot == 0) return ray_plane(T, D, N, d);

      //Create the rest of the variables
      double rdot = T | N;
      
      //Ensure that the normal is pointing towards the particle
      //position (so that (rdot < d) is a test if the particle is
      //overlapped)
      if (rdot < 0)
	{ adot = -adot; rdot = -rdot; N = -N; }
      
      double vdot = D | N;

      //Check for overlapped and approaching dynamics
      if ((rdot < d) && (vdot < 0)) return 0;

      double arg = vdot * vdot - 2 * (rdot - d) * adot;
      double minimum = - vdot / adot;

      if (adot < 0)
	{ //Particle will always hit the wall
	  
	  //If there are no real roots, the particle is arcing inside
	  //the wall and there is a collision at the minimum
	  if (arg < 0) return minimum;
	  
	  std::pair<double, double> roots
	    = magnet::math::quadraticEquation(0.5 * adot, vdot, rdot - d);
	  return std::max(roots.first, roots.second);
	}
      else
	{ //Particle can curve away from the wall

	  //Check if the particle misses the wall completely or we're
	  //passed the minimum
	  if ((arg < 0) || (minimum < 0)) return HUGE_VAL;
	  
	  std::pair<double, double> roots
	    = magnet::math::quadraticEquation(0.5 * adot, vdot, rdot - d);
	  return std::min(roots.first, roots.second);	  
	}
    }
  }
}
