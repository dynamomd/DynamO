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
#include "trajectory.hpp"
#include "../../dynamics/dynamics.hpp"
#include "../../dynamics/units/units.hpp"
#include "../../dynamics/globals/globEvent.hpp"
#include "../../dynamics/interactions/intEvent.hpp"
#include "../../dynamics/locals/localEvent.hpp"
#include "../../dynamics/NparticleEventData.hpp"
#include "../../dynamics/systems/system.hpp"
#include "../../dynamics/BC/BC.hpp"

COPTrajectory::COPTrajectory(const DYNAMO::SimData* t1, const XMLNode&):
  COutputPlugin(t1,"Trajectory")
{}

COPTrajectory::COPTrajectory(const COPTrajectory& trj):
  COutputPlugin(trj)
{
  if (trj.logfile.is_open())
    trj.logfile.close();
}

void
COPTrajectory::initialise()
{
  if (logfile.is_open())
    logfile.close();
  
  logfile.open("trajectory.out", std::ios::out|std::ios::trunc);
}

void
COPTrajectory::printData(const CParticle& p1,
			 const CParticle& p2) const
{
  size_t id1 = ((p1.getID() < p2.getID()) 
		? p1.getID() : p2.getID());
  
  size_t id2 = ((p1.getID() > p2.getID()) 
		? p1.getID() : p2.getID());

  CVector<> rij = Sim->vParticleList[id1].getPosition()
    - Sim->vParticleList[id2].getPosition();

  Sim->Dynamics.BCs().setPBC(rij);
  
  rij /= Sim->Dynamics.units().unitLength();

  logfile << " p1 " << id1
	  << " p2 " << id2
	  << " |r12| < " << rij.length()
	  << " r12 < ";
  
  for (size_t iDim(0); iDim < NDIM; ++iDim)
    logfile << rij[iDim] << " ";
}

void 
COPTrajectory::eventUpdate(const CIntEvent& eevent, 
				   const C2ParticleData&)
{
  logfile << "INTERACTION " << eevent.getInteraction().getID()
	  << " t " << Sim->dSysTime / Sim->Dynamics.units().unitTime() 
	  << " dt " << eevent.getdt() / Sim->Dynamics.units().unitTime();
  
  printData(eevent.getParticle1(),
	    eevent.getParticle2());
  
  logfile << "\n";
}

void 
COPTrajectory::eventUpdate(const CGlobEvent& eevent, 
				   const CNParticleData& SDat)
{
  logfile << "GLOBAL " << eevent.getGlobalID()
	  << " t " << Sim->dSysTime / Sim->Dynamics.units().unitTime() 
	  << " dt " << eevent.getdt() / Sim->Dynamics.units().unitTime()
	  << "\n";

  BOOST_FOREACH(const C1ParticleData& pData, SDat.L1partChanges)
    logfile << "    1PEvent p1" << pData.getParticle().getID()
	    << "\n";
  
  BOOST_FOREACH(const C2ParticleData& pData, SDat.L2partChanges)
    {
      logfile << "    2PEvent";

      printData(pData.particle1_.getParticle(),
		pData.particle2_.getParticle());

      logfile << "\n";
    }
}

void 
COPTrajectory::eventUpdate(const CLocalEvent& eevent, 
				   const CNParticleData& SDat)
{
  logfile << "LOCAL " << eevent.getLocalID()
	  << " t " << Sim->dSysTime / Sim->Dynamics.units().unitTime() 
	  << " dt " << eevent.getdt() / Sim->Dynamics.units().unitTime()
	  << "\n";

  BOOST_FOREACH(const C1ParticleData& pData, SDat.L1partChanges)
    logfile << "    1PEvent p1 " << pData.getParticle().getID()
	    << "\n";
  
  BOOST_FOREACH(const C2ParticleData& pData, SDat.L2partChanges)
    {
      logfile << "    2PEvent";
      
      printData(pData.particle1_.getParticle(),
		pData.particle2_.getParticle());
      
      logfile << "\n";
    }
}

void 
COPTrajectory::eventUpdate(const CSystem& sys, const CNParticleData& SDat, 
				   const Iflt& dt)
{
  logfile << "SYSTEM " << sys.getID()
	  << " t " << Sim->dSysTime / Sim->Dynamics.units().unitTime() 
	  << " dt " << dt / Sim->Dynamics.units().unitTime()
	  << "\n";

  BOOST_FOREACH(const C1ParticleData& pData, SDat.L1partChanges)
    logfile << "    1PEvent p1 " << pData.getParticle().getID()
	    << "\n";
  
  BOOST_FOREACH(const C2ParticleData& pData, SDat.L2partChanges)
    {
      logfile << "    2PEvent";
      
      printData(pData.particle1_.getParticle(),
		pData.particle2_.getParticle());
      
      logfile << "\n";
    }
}

void 
COPTrajectory::output(xmlw::XmlStream& XML)
{}
