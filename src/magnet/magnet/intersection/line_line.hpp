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
#include <magnet/intersection/overlapfuncs/lines.hpp>
#include <magnet/math/frenkelroot.hpp>


namespace magnet {
  namespace intersection {
    /*! \brief A line-line intersection test.
     */
    inline std::pair<bool, double> line_line(const math::Vector& rij, const math::Vector& vij,
					     const math::Vector& angvi, const math::Vector& angvj,
					     const math::Quaternion& orientationi, const math::Quaternion& orientationj,
					     const double& length, bool skip_zero, double t_max)
    {
      overlapfuncs::Lines fL(rij, vij, angvi, angvj, orientationi, orientationj, length);
      
      //Shift the lower bound up so we don't find the same root again
      double t_min = skip_zero ? fabs(2.0 * fL.eval<1>()) / fL.max<2>() : 0;
    
      //Find window delimited by discs
      std::pair<double,double> dtw = fL.discIntersectionWindow();
      t_min = std::max(dtw.first, t_min);
      t_max = std::min(dtw.second, t_max);
      return magnet::math::frenkelRootSearch(fL, t_min, t_max, length * 1e-10);
    }
  }
}
