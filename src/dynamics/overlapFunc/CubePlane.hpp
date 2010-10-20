/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "../../base/constants.hpp"
#include "../../datatypes/vector.hpp"

#include <limits>

namespace DYNAMO
{
  namespace OverlapFunctions
  {
    /*! \brief Discovers if a cube and a plane intersect by testing
     * which side of the plane the points of the cube lie on.
     *
     * This is used in the collision CLSentinel to install itself in cells
     */
    bool CubePlane(const Vector& CubeOrigin, const Vector& CubeDimensions,
		   const Vector& PlaneOrigin, const Vector& PlaneNormal,
		   const double tol = 0);
  };
};
