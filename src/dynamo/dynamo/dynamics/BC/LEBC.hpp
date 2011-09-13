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
#include <dynamo/dynamics/BC/BC.hpp>

namespace dynamo {
  class Particle;

  /*! \brief A Lees-Edwards simple shear boundary condition.
   * 
   * This class implements the sliding brick boundary condtion. In this
   * the simulation images above and below the primary image are set in
   * motion. This affects the particle velocities and positions on a
   * transition of the boundary.
   *
   * See \ref BoundaryCondition for a general description of the member functions.
   */
  class BCLeesEdwards: public BoundaryCondition
  {
  public:
    BCLeesEdwards(const dynamo::SimData*);

    BCLeesEdwards(const magnet::xml::Node&, const dynamo::SimData*);

    virtual void outputXML(magnet::xml::XmlStream&) const;

    virtual void operator<<(const magnet::xml::Node&);

    virtual void applyBC(Vector&) const; 

    virtual void applyBC(Vector&, Vector&) const;

    virtual void applyBC(Vector&, const double& dt) const;

    virtual void update(const double&);

    /*! \brief Returns the shear rate of the boundaries. */
    inline double getShearRate() const { return _shearRate; }

    /*! \brief Returns the stream velocity at the passed particles
     *    position. 
     *
     * This stream velocity is based off a linear interpolation between
     * the boundary velocities. It is only guaranteed to be correct at
     * the simulation boundaries and its periodic images.
     *
     * This is an important distinction to make, as Lees-Edwards
     * boundary conditions do not force a linear shear profile, only a
     * specified fixed shear rate over a box length. To enforce a linear
     * profile you must also use a thermostat, but this may be
     * problematic (see Evans and Morris, "Statistical Mechanics of
     * Nonequilibrium Liquids"). Essentially, a thermostat will cause
     * "strings" to form in the system.
     */
    Vector getStreamVelocity(const Particle& part) const;

    /*! \brief Returns the peculiar velocity of the particle.
     *
     * By definition, the peculiar velocity is the velocity of a
     * particle, minus the velocity of the fluid at that point.  \sa
     * getStreamVelocity(const Particle&)
     */
    Vector getPeculiarVelocity(const Particle& part) const;

  protected:  
    /*! \brief The amount neighboring periodic images have slid against
     *   each other.
     * 
     * This value must be stored so that when a simulation is saved and
     * loaded so the sliding PBC images are at the same place.
     */
    double _dxd;

    /*! \brief The rate of shear.*/
    double _shearRate;
  };
}
