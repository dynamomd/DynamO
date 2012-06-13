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
#include <dynamo/globals/global.hpp>
#include <vector>

namespace dynamo {
  /*! \brief A Global event which helps prevent wrap around problems
   * with neighbor lists in periodic systems.
   *
   * If a particle has a clear path from one end of the simulation to
   * the other and PBC are applied, the cellular neighbor lists can
   * enter an infinite loop. The particle keeps travelling around and
   * around the simulation, without actually moving forward in time
   * because it doesn't actually hit anything.
   *
   * This Global attempts to break this infinite loop by making
   * particles have a virtual event when they travel half a simulation
   * box length.
   */
  class GPBCSentinel: public Global
  {
  public:
    GPBCSentinel(const magnet::xml::Node&, dynamo::Simulation*);

    GPBCSentinel(dynamo::Simulation*, const std::string&);
  
    virtual ~GPBCSentinel() {}

    virtual GlobalEvent getEvent(const Particle &) const;

    virtual void runEvent(Particle&, const double) const;

    virtual void initialise(size_t);

    virtual void operator<<(const magnet::xml::Node&);

  protected:
    virtual void outputXML(magnet::xml::XmlStream&) const {}

    double maxintdist;
  };
}
