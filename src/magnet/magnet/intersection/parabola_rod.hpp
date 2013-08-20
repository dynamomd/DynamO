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
#include <magnet/intersection/parabola_cylinder.hpp>

namespace magnet {
  namespace intersection {
    /*! \brief A parabola-rod intersection test.
      
      A rod is a cylinder which is not infinite, but of limited
      length. The cylinder is defined using a single base vertex at the
      center of the bottom circular face and an axial vector pointing
      from the base vertex to the top vertex. This test ignores the
      back face of the rod. It is used to detect when a parabolic ray
      will enter a rod.
     
      \param T The origin of the ray relative to the base vertex.
      \param D The direction/velocity of the ray.
      \param Aray The acceleration acting on the ray.
      \param A The axial vector of the rod.
      \param r Radius of the rod.
      \return The time until the intersection, or HUGE_VAL if no intersection.
    */
    inline double parabola_rod(math::Vector T, math::Vector D, math::Vector Aray, const math::Vector& A, const double r)
    {
      double t = parabola_cylinder(T, D, Aray, A / A.nrm(), r);
      double Tproj = ((T + t * D + 0.5 * t * t * Aray) | A);
      if ((Tproj < 0) || (Tproj > A.nrm2())) return HUGE_VAL;
      return t;
    }
  }
}
