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

#ifndef CPDData_H
#define CPDData_H
#include "../../datatypes/vector.hpp"
#include "../../base/is_simdata.hpp"
#include "../dynamics.hpp"
#include "../BC/BC.hpp"

//Pair dynamics data, can be allocated in compile time saving new() calls
struct CPDData
{
  inline CPDData() {}

  inline CPDData(const DYNAMO::SimData& Sim, const CParticle& p1, 
		 const CParticle& p2):
    rij(p1.getPosition() - p2.getPosition()),
    vij(p1.getVelocity() - p2.getVelocity()),    
    dt(HUGE_VAL)
  {
    Sim.Dynamics.BCs().applyBC(rij, vij);
    rvdot = rij | vij;
    r2 = rij.nrm2();
    v2 = vij.nrm2();
  }

  CPDData(const DYNAMO::SimData&, const CRange&, const CRange&);

  Vector  rij, vij;
  Iflt rvdot;
  Iflt r2;
  Iflt v2;
  Iflt dt;
};

#endif
