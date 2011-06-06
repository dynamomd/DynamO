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
#include "../../datatypes/vector.hpp"
#include "../../base/is_simdata.hpp"
#include "../dynamics.hpp"
#include "../BC/BC.hpp"

//Pair dynamics data, can be allocated in compile time saving new() calls
struct CPDData
{
  inline CPDData(): p1(NULL), p2(NULL) {}

  inline CPDData(const dynamo::SimData& Sim, const Particle& np1, 
		 const Particle& np2):
    rij(np1.getPosition() - np2.getPosition()),
    vij(np1.getVelocity() - np2.getVelocity()),    
    dt(HUGE_VAL),
    p1(&np1),
    p2(&np2)
  {
    Sim.dynamics.BCs().applyBC(rij, vij);
    rvdot = rij | vij;
    r2 = rij.nrm2();
    v2 = vij.nrm2();
  }

  CPDData(const dynamo::SimData&, const CRange&, const CRange&);

  Vector  rij, vij;
  double rvdot;
  double r2;
  double v2;
  double dt;
  const Particle* const p1;
  const Particle* const p2;
};
