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
#include <dynamo/dynamics/interactions/interaction.hpp>
#include <magnet/math/vector.hpp>

namespace dynamo {
  //! This class is inherited by Interaction -s that can be represented
  //! as a collection of spheres.
  //!
  //! The methods provided in this base class allow a spherical representation to be created by a Species to pass to the coil visualizer. 
  class SphericalRepresentation
  {
  public:

    //! This function returns how many spheres must be rendered per
    //! particle. Typically it is only one but Interaction classes like
    //! IDumbbells need two or more.
    virtual size_t spheresPerParticle() const = 0;

    //! Returns the diameter of the one of the spheres used to represent
    //! a particle.
    //! \param ID The id of the particle to fetch the diameter for.
    //! \param subID The index of one of the spheres used to render the particle. This must be in the range [0, spheresPerParticle()-1].
    virtual double getDiameter(size_t ID, size_t subID) const = 0;

    //! Returns the location of the one of the spheres used to represent
    //! a particle.
    //! \param ID The id of the particle to fetch the location for.
    //! \param subID The index of one of the spheres used to render the particle. This must be in the range [0, spheresPerParticle()-1].
    virtual Vector getPosition(size_t ID, size_t subID) const = 0;
  };
}

