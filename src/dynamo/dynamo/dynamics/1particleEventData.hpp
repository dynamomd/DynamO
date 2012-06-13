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
#include <dynamo/dynamics/species/species.hpp>
#include <dynamo/dynamics/eventtypes.hpp>
#include <dynamo/simulation/particle.hpp>

namespace dynamo {
  class ParticleEventData
  {
  public:
    ParticleEventData(const Particle& part, const Species& sp, 
		      EEventType eType):
      particle_(part), oldVelVec(part.getVelocity()),
      species_(sp), Type_(eType), 
      deltaU(0.0), deltaKE(0.0)
    {}

    inline const Particle& getParticle() const
    { return particle_; }
  
    inline const Vector & getOldVel() const
    { return oldVelVec; }

    inline Vector  getOldPosition() const
    { M_throw() << "Not yet Implemented"; }

    inline const Species& getSpecies() const
    { return species_; }

    inline void setType(EEventType nType)
    { Type_ = nType; }

    inline double getDeltaU() const
    { return deltaU; }
  
    inline void setDeltaU(const double& dU)
    { deltaU = dU; }

    inline double getDeltaKE() const
    { return deltaKE; }
  
    inline void setDeltaKE(const double& dKE)
    { deltaKE = dKE; }

    inline EEventType getType() const
    { return Type_; }
  
    Vector  getDeltaP() const
    { return species_.getMass(particle_.getID()) * (particle_.getVelocity() - oldVelVec); }

  private:
    const Particle& particle_;
    const Vector  oldVelVec;
    const Species& species_;
    EEventType Type_;
    double deltaU;
    double deltaKE;
  };
}
