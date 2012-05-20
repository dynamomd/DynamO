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

#include <dynamo/dynamics/1particleEventData.hpp>
#include <dynamo/simulation/particle.hpp>

namespace dynamo {
  class PairEventData
  {
  public:
    PairEventData(const Particle& part1,
		  const Particle& part2,
		  const Species& sp1,
		  const Species& sp2,
		  EEventType eType):
      particle1_(part1, sp1, eType), 
      particle2_(part2, sp2, eType),
      rij(part1.getPosition() - part2.getPosition()),
      vijold(part1.getVelocity() - part2.getVelocity())
    {}
  
    ParticleEventData particle1_;
    ParticleEventData particle2_;
    Vector  rij;
    Vector  vijold;
    Vector  dP;
    double rvdot;

    void setType(EEventType nType)
    { particle1_.setType(nType); particle2_.setType(nType); }

    EEventType getType() const
    { return particle1_.getType(); }
  };
}
