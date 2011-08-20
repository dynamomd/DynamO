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

#include <dynamo/dynamics/liouvillean/datastruct.hpp>
#include <dynamo/dynamics/include.hpp>
#include <boost/foreach.hpp>

namespace dynamo {
  CPDData::CPDData(const dynamo::SimData& Sim, const CRange& range1, 
		   const CRange& range2):
    dt(HUGE_VAL),
    p1(NULL),
    p2(NULL)
  {
    Vector  COMVel1(0,0,0), COMVel2(0,0,0), COMPos1(0,0,0), COMPos2(0,0,0);
  
    double structmass1(0), structmass2(0);
    
    BOOST_FOREACH(const size_t& ID, range1)
      {
	double mass = Sim.dynamics.getSpecies(Sim.particleList[ID]).getMass(ID);

	structmass1 += mass;
	COMVel1 += Sim.particleList[ID].getVelocity() * mass;
	COMPos1 += Sim.particleList[ID].getPosition() * mass;
      }
    
    BOOST_FOREACH(const size_t& ID, range2)
      {
	double mass = Sim.dynamics.getSpecies(Sim.particleList[ID]).getMass(ID);
	structmass2 += mass;
	COMVel2 += Sim.particleList[ID].getVelocity() * mass;	
	COMPos2 += Sim.particleList[ID].getPosition() * mass;
      }
    
    COMVel1 /= structmass1;
    COMVel2 /= structmass2;

    COMPos1 /= structmass1;
    COMPos2 /= structmass2;

    rij = COMPos1 - COMPos2;

    vij = COMVel1 - COMVel2;

    Sim.dynamics.BCs().applyBC(rij, vij);

    rvdot = (rij | vij);

    r2 = rij.nrm2();
    v2 = vij.nrm2();
  }
}
