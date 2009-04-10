/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef C1ParticleData_HPP
#define C1ParticleData_HPP

#include "../datatypes/vector.hpp"
#include "species/species.hpp"
#include "../simulation/particle.hpp"

class C1ParticleData
{
public:
  C1ParticleData(const CParticle& part, const CSpecies& sp, 
		 EEventType eType):
    particle_(part), oldVelVec(part.getVelocity()),
    species_(sp), Type_(eType), 
    deltaU(0.0), deltaKE(0.0)
  {}

  inline const CParticle& getParticle() const
  { return particle_; }
  
  inline const Vector & getOldVel() const
  { return oldVelVec; }

  inline Vector  getOldPosition() const
  { D_throw() << "Not yet Implemented"; }

  inline const CSpecies& getSpecies() const
  { return species_; }

  inline void setType(EEventType nType)
  { Type_ = nType; }

  inline Iflt getDeltaU() const
  { return deltaU; }
  
  inline void setDeltaU(const Iflt& dU)
  { deltaU = dU; }

  inline Iflt getDeltaKE() const
  { return deltaKE; }
  
  inline void setDeltaKE(const Iflt& dKE)
  { deltaKE = dKE; }

  inline EEventType getType() const
  { return Type_; }
  
  Vector  getDeltaP() const
  {
    return species_.getMass() * (particle_.getVelocity() - oldVelVec);
  }

private:
  const CParticle& particle_;
  const Vector  oldVelVec;
  const CSpecies& species_;
  EEventType Type_;
  Iflt deltaU;
  Iflt deltaKE;
};

#endif
