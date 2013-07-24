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
#include <dynamo/species/species.hpp>
#include <dynamo/eventtypes.hpp>
#include <dynamo/particle.hpp>

namespace dynamo {
  class ParticleEventData
  {
  public:
    ParticleEventData() {}

    ParticleEventData(const Particle& part, const Species& sp,
		      EEventType eType):
      particleID(part.getID()), 
      oldVelVec(part.getVelocity()),
      speciesID(sp.getID()), Type_(eType), 
      deltaU(0.0)
    {}

    inline size_t getParticleID() const
    { return particleID; }
  
    inline Vector getOldVel() const
    { return oldVelVec; }

    inline Vector getOldPosition() const
    { M_throw() << "Not yet Implemented"; }

    inline size_t getSpeciesID() const
    { return speciesID; }

    inline void setType(EEventType nType)
    { Type_ = nType; }

    inline double getDeltaU() const
    { return deltaU; }
  
    inline void setDeltaU(const double& dU)
    { deltaU = dU; }

    inline EEventType getType() const
    { return Type_; }
  
  private:
    size_t particleID;
    Vector  oldVelVec;
    size_t speciesID;
    EEventType Type_;
    double deltaU;
  };
}
