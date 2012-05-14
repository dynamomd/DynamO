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
    /*! \brief A ray->Axis-Aligned-Cube intersection test with back
        face culling.

      This test works by determining when the ray entered the axis
      aligned "range" of each cubes face. Then it is tested if the
      entries are all before the exits, if true the ray intersects the
      cube at some point. Finally, we check if the intersecting ray is
      approaching the cube at the entry point in order to cull back
      faces.
      
      \param T The origin of the ray relative to the cube center.
      \param D The direction/velocity of the ray.
      \param C The dimensions of the cube.
      \return The time until the intersection, or HUGE_VAL if no intersection.
    */
    inline double ray_AAcube_bfc(const math::Vector& T,
				 const math::Vector& D,
				 math::Vector C)
    {
      //We need the cube half lengths
      C *= 0.5;

      double time_in_max = -HUGE_VAL;
      double time_out_min = HUGE_VAL;
      
      for (size_t i(0); i< 3; ++i)
	{
	  //Test if the velocity is zero
	  if (D[i] == 0)
	    {//Its zero, a collision can only occur if this range is already within the cubes limits.
	      if ((T[i] * T[i]) > (C[i] * C[i]))
		return HUGE_VAL;
	      else
		{
		  //We're within the cube limits but with zero
		  //velocity, this means we entered at -HUGE_VAL and
		  //left at HUGE_VAL. Nothing to be done here
		}
	    }
	  else
	    {
	      double time_in  = (-copysign(C[i],D[i]) - T[i]) / D[i];
	      double time_out = (+copysign(C[i],D[i]) - T[i]) / D[i];

	      time_in_max = std::max(time_in_max, time_in);
	      time_out_min = std::min(time_out_min, time_out);
	    }
	}
      
      //Now we check if an intersection occurs
      if (time_in_max > time_out_min)
	return HUGE_VAL;
      
      //An intersection occurs, we only cause an intersection if the
      //ray is still entering the object (the ray is closer to the entry
      //than the exit).
      if (std::abs(time_in_max) < std::abs(time_out_min))
	return time_in_max;
      
      //There is an intersection, but the ray is already leaving the
      //cube.
      return HUGE_VAL;
    }
  }
}
