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
  namespace overlap {
    //! \brief A point-sphere overlap test.
    //!
    //! This function assumes the point location passed is relative to
    //! the sphere's center.
    //!
    //! \param S The point's relative position.
    //! \param d The diameter of the sphere.
    //! \return Whether the point is inside the sphere.
    inline bool point_sphere(const math::Vector& P, 
			     const double d)
    { return P.nrm2() <= d * d; }

    //! \brief A sphere-sphere overlap test.
    //!
    //! This function assumes the sphere location passed is relative to
    //! the other sphere's center.
    //!
    //! \param S The other sphere's relative position.
    //! \param d The average diameter of the spheres.
    //! \return Whether the spheres are overlapping.
    inline bool sphere_sphere(const math::Vector& P, 
			     const double d)
    { return point_sphere(P, d); }
  }
}
