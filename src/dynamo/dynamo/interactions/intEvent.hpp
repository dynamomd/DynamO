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

#include <dynamo/eventtypes.hpp>
#include <dynamo/particle.hpp>
#include <dynamo/interactions/interaction.hpp>
#include <cfloat>
#include <limits>

struct XMLNode;

class Interaction;
namespace dynamo {
  class Simulation;
}

namespace xml
{
  class XmlStream;
}

namespace dynamo {
  class IntEvent
  {
  public:  
    friend class Event;

    inline IntEvent ():
      particle1(std::numeric_limits<size_t>::max()), particle2(std::numeric_limits<size_t>::max()), 
      dt(HUGE_VAL), CType(NONE),
      intID(std::numeric_limits<size_t>::max()) {}

    inline IntEvent(const Particle& part1, const Particle& part2, 
		    const double &delt, EEventType nType, 
		    const Interaction& pI):
      particle1(part1.getID()), particle2(part2.getID()), dt(delt), 
      CType(nType), intID(pI.getID()) {}
  
    inline IntEvent (const Particle& part1):
      particle1(part1.getID()), particle2(std::numeric_limits<size_t>::max()), 
      dt(HUGE_VAL), CType(NONE),
      intID(std::numeric_limits<size_t>::max()) {}

    inline IntEvent(const Particle& part1, const double& dt, 
		    EEventType etype):
      particle1(part1.getID()), particle2(std::numeric_limits<size_t>::max()), 
      dt(dt), CType(etype),
      intID(std::numeric_limits<size_t>::max()) {}

    inline bool operator== (const Particle &partx) const 
    { return ((particle1 == partx.getID()) || (particle2 == partx.getID())); }
  
    inline bool areInvolved(const IntEvent &coll) const 
    { 
      return ((coll.particle1 == particle1) 
	      || (coll.particle1 == particle2)
	      || (coll.particle2 == particle1)
	      || (coll.particle2 == particle2));
    }
  
    inline void invalidate() 
    { 
      dt = DBL_MAX; 
      CType = NONE; 
    }

    inline bool operator< (const IntEvent & C2) const { return dt < C2.dt;}
    inline bool operator> (const IntEvent & C2) const { return dt > C2.dt;}
    inline void incrementTime(const double deltat) {dt -= deltat; }
    inline void addTime(const double deltat) {dt += deltat; }
    inline const size_t& getParticle1ID() const { return particle1; }
    inline const size_t& getParticle2ID() const { return particle2; }
    inline bool hasParticle2() const { return particle2 != std::numeric_limits<size_t>::max(); }
    inline const double& getdt() const { return dt; }
    inline EEventType getType() const { return CType; }
    inline void setType(EEventType a) const { CType = a; }
    inline void scaleTime(const double& scale) { dt *= scale; }
    inline const size_t& getInteractionID() const { return intID; }

  private:
    size_t  particle1;
    size_t  particle2;
    double dt;
    mutable EEventType CType;
    size_t intID;
  };
}
