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

#ifndef C2ParticleData_HPP
#define C2ParticleData_HPP

#include "1particleEventData.hpp"
#include "../datatypes/vector.hpp"
#include "../simulation/particle.hpp"

class C2ParticleData
{
public:
  C2ParticleData(const CParticle& part1,
		 const CParticle& part2,
		 const CSpecies& sp1,
		 const CSpecies& sp2,
		 EEventType eType):
    particle1_(part1, sp1, eType), 
    particle2_(part2, sp2, eType),
    rij(part1.getPosition() - part2.getPosition()),
    vijold(part1.getVelocity() - part2.getVelocity())
  {}
  
  C1ParticleData particle1_;
  C1ParticleData particle2_;
  CVector<> rij;
  CVector<> vijold;
  CVector<> dP;
  Iflt rvdot;

  void setType(EEventType nType)
  { particle1_.setType(nType); particle2_.setType(nType); }

  EEventType getType() const
  { return particle1_.getType(); }
};

#endif
