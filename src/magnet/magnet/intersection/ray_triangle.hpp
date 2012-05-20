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
    //! \brief A ray-triangle intersection test.
    //!
    //! The method used is described in "Fast, minimum storage
    //! ray-triangle intersection", by Tomas Möller and Ben Trumbore.
    //!
    //! Here we make the assumption that all positions are relative to
    //! the first vertex. Thus we only need to pass the two edge
    //! vectors of the triangle.
    //!
    //! For backface culling, we use counter-clockwise ordering of the
    //! vertices.
    //!
    //! \tparam BACKFACE_CULLING Ignores ray triangle intersections
    //! where the ray enters the back face of the triangle.
    //! \tparam DIAGONAL_TEST Enables one of the checks used to make
    //! sure the intersection point is within the triangle. This
    //! should only be turned off to implement other ray-shape
    //! intersection tests (\sa ray_quadrilateral).
    //!
    //! \param T The origin of the ray relative to the first vertex.
    //! \param D The direction/velocity of the ray.
    //! \param E1 The first edge vector of the triangle (V1-V0).
    //! \param E2 The second edge vector of the triangle (V2-V0).
    //! \return The time until the intersection, or HUGE_VAL if no intersection.
    template<bool BACKFACE_CULLING, bool DIAGONAL_TEST>
    inline double ray_triangle(const math::Vector& T, 
			       const math::Vector& D,
			       const math::Vector& E1, 
			       const math::Vector& E2)
    {
      math::Vector P = D ^ E2;
      double det = E1 | P;
      
      if (BACKFACE_CULLING)
	{
	  //Ray is parallel (0) or the ray is leaving the triangle (not entering it)
	  if (det <= 0) return HUGE_VAL;
	  
	  double u = T | P;
	  if ((u < 0) || (u > det)) return HUGE_VAL;
	  
	  math::Vector Q = T ^ E1;
	  double v = D | Q;
	  
	  if ((v < 0) 
	      || (DIAGONAL_TEST && ((u + v) > det)))  return HUGE_VAL;
	  
	  return (E2 | Q) / det;
	}
      else
	{
	  if (det == 0) return HUGE_VAL;
	  
	  //We must work out the inv determinate now and use it in the
	  //calculations as we don't know it's sign.
	  double invdet = 1/ det;

	  double u = (T | P) * invdet;

	  if ((u < 0) || (u > 1)) return HUGE_VAL;

	  math::Vector Q = T ^ E1;
	  double v = (D | Q) * invdet;

	  if ((v < 0) || (DIAGONAL_TEST && ((u + v) > 1)))  return HUGE_VAL;

	  return (E2 | Q) * invdet;
	}
    }
  }
}
