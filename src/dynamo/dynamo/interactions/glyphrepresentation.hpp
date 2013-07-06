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
#include <dynamo/interactions/interaction.hpp>
#include <magnet/math/vector.hpp>

namespace dynamo {
  /*! \brief This class is inherited by Interaction -s that can be
    represented as a collection of spheres.
  
    The methods provided in this base class allow a spherical
    representation to be created by a Species to pass to the coil
    visualizer.
  */
  class GlyphRepresentation { 
  public:
    enum GLYPH_TYPE
      {
	SPHERE_GLYPH=0,
	ARROW_GLYPH=1,
	CYLINDER_GLYPH=2,
	LINE_GLYPH=3,
	CUBE_GLYPH=4
      };

    /*! Returns the diameter of the one of the spheres used to represent
      a particle.
      
      \param ID The id of the particle to fetch the diameter for.
      
      \param subID The index of one of the spheres used to render the
      particle. This must be in the range [0, spheresPerParticle()-1].
    */
    virtual Vector getGlyphSize(size_t ID) const = 0;

    /*! Returns the location of the one of the spheres used to
      represent a particle.
      
      \param ID The id of the particle to fetch the location for.
      
      \param subID The index of one of the spheres used to render the
      particle. This must be in the range [0, spheresPerParticle()-1].
    */
    virtual GLYPH_TYPE getDefaultGlyphType() const { return SPHERE_GLYPH; }
  };
}

