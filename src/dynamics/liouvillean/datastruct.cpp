#include "datastruct.hpp"
#include <boost/foreach.hpp>
#include "../include.hpp"

CPDData::CPDData(const DYNAMO::SimData& Sim, const CRange& range1, 
		 const CRange& range2):
  dt(HUGE_VAL),
  p1(NULL),
  p2(NULL)
{
  Vector  COMVel1(0,0,0), COMVel2(0,0,0), COMPos1(0,0,0), COMPos2(0,0,0);
  
  double structmass1(0), structmass2(0);
    
  BOOST_FOREACH(const size_t& ID, range1)
    {
      structmass1 += 
	Sim.dynamics.getSpecies(Sim.particleList[ID]).getMass();
	
      COMVel1 += Sim.particleList[ID].getVelocity()
	* Sim.dynamics.getSpecies(Sim.particleList[ID]).getMass();
	
      COMPos1 += Sim.particleList[ID].getPosition()
	* Sim.dynamics.getSpecies(Sim.particleList[ID]).getMass();
    }
    
  BOOST_FOREACH(const size_t& ID, range2)
    {
      structmass2 += 
	Sim.dynamics.getSpecies(Sim.particleList[ID]).getMass();
	
      COMVel2 += Sim.particleList[ID].getVelocity()
	* Sim.dynamics.getSpecies(Sim.particleList[ID]).getMass();
	
      COMPos2 += Sim.particleList[ID].getPosition()
	* Sim.dynamics.getSpecies(Sim.particleList[ID]).getMass();
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
