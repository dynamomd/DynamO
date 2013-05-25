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
#include <dynamo/outputplugins/general/trajectory.hpp>

#include <dynamo/units/units.hpp>
#include <dynamo/globals/globEvent.hpp>
#include <dynamo/interactions/intEvent.hpp>
#include <dynamo/locals/localEvent.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/systems/system.hpp>
#include <dynamo/BC/BC.hpp>
#include <iomanip>

namespace dynamo {
  OPTrajectory::OPTrajectory(const dynamo::Simulation* t1, const magnet::xml::Node&):
    OutputPlugin(t1,"Trajectory")
  {}

  OPTrajectory::OPTrajectory(const OPTrajectory& trj):
    OutputPlugin(trj)
  {
    if (trj.logfile.is_open())
      trj.logfile.close();
  }

  void
  OPTrajectory::initialise()
  {
    if (logfile.is_open())
      logfile.close();
 
    logfile.open("trajectory.out", std::ios::out|std::ios::trunc);
    logfile.precision(4);
    logfile.setf(std::ios::fixed);
  }

  void
  OPTrajectory::printData(const size_t& p1,
			  const size_t& p2) const
  {
    size_t id1 = ((p1 < p2) 
		  ? p1 : p2);
  
    size_t id2 = ((p1 > p2) 
		  ? p1 : p2);

    Vector  rij = Sim->particles[id1].getPosition()
      - Sim->particles[id2].getPosition(),
      vij = Sim->particles[id1].getVelocity()
      - Sim->particles[id2].getVelocity();
  

    Sim->BCs->applyBC(rij, vij);
  
    rij /= Sim->units.unitLength();
    vij /= Sim->units.unitVelocity();

    logfile << " p1 " << std::setw(5) << id1
	    << " p2 " << std::setw(5) << id2
	    << " |r12| " << std::setw(5) << rij.nrm()
	    << " post-r12 < ";
  
    for (size_t iDim(0); iDim < NDIM; ++iDim)
      logfile << std::setw(7) << rij[iDim] << " ";

    logfile << ">";

    logfile << " post-v12 < ";
    for (size_t iDim(0); iDim < NDIM; ++iDim)
      logfile << std::setw(7) << vij[iDim] << " ";

    logfile << "> post-rvdot " << (vij | rij) ;
  }

  void 
  OPTrajectory::eventUpdate(const IntEvent& eevent, 
			    const PairEventData& pdat)
  {
    logfile << std::setw(8) << Sim->eventCount
	    << " INTERACTION " << eevent.getInteractionID()
	    << " TYPE " << eevent.getType()
	    << " t " << std::setw(5) << Sim->systemTime / Sim->units.unitTime() 
	    << " dt " << std::setw(5) << eevent.getdt() / Sim->units.unitTime();

    logfile << " deltaP1 < ";
    for (size_t iDim(0); iDim < NDIM; ++iDim)
      logfile << std::setw(7) 
	      << ((eevent.getParticle1ID() < eevent.getParticle2ID())? -1 : 1) 
	* pdat.impulse[iDim] << " ";

    logfile << " >";
  
    printData(eevent.getParticle1ID(),
	      eevent.getParticle2ID());
  
    logfile << "\n";
  }

  void 
  OPTrajectory::eventUpdate(const GlobalEvent& eevent, 
			    const NEventData& SDat)
  {
    logfile << std::setw(8) << Sim->eventCount
	    << " GLOBAL " << eevent.getGlobalID()
	    << " TYPE " << eevent.getType()
	    << " t " << Sim->systemTime / Sim->units.unitTime() 
	    << " dt " << eevent.getdt() / Sim->units.unitTime()
	    << "\n";

    for (const ParticleEventData& pData : SDat.L1partChanges)
      {
	const Particle& part = Sim->particles[pData.getParticleID()];
	logfile << "    1PEvent p1 " << part.getID();
	Vector delP = Sim->species[pData.getSpeciesID()]->getMass(part.getID()) * (part.getVelocity() - pData.getOldVel());
	delP /= Sim->units.unitMomentum();
	Vector pos = part.getPosition() / Sim->units.unitLength();
	Vector oldv = pData.getOldVel() / Sim->units.unitVelocity();
	Vector newv = part.getVelocity() / Sim->units.unitVelocity();
	logfile << " delP1=" << delP.toString() << ", pos=" << pos.toString() << ", vel=" << newv.toString() << ", oldvel=" << oldv.toString() << "\n";
      }
  
    for (const PairEventData& pData : SDat.L2partChanges)
      {
	logfile << "    2PEvent";

	printData(pData.particle1_.getParticleID(),
		  pData.particle2_.getParticleID());

	logfile << "\n";
      }
  }

  void 
  OPTrajectory::eventUpdate(const LocalEvent& eevent, 
			    const NEventData& SDat)
  {
    logfile << std::setw(8) << Sim->eventCount 
	    << " LOCAL " << eevent.getLocalID()
	    << " TYPE " << eevent.getType()
	    << " t " << Sim->systemTime / Sim->units.unitTime() 
	    << " dt " << eevent.getdt() / Sim->units.unitTime()
	    << "\n";

    for (const ParticleEventData& pData : SDat.L1partChanges)
      {
	const Particle& part = Sim->particles[pData.getParticleID()];
	logfile << "    1PEvent p1 " << part.getID();
	Vector delP = Sim->species[pData.getSpeciesID()]->getMass(part.getID()) * (part.getVelocity() - pData.getOldVel());
	delP /= Sim->units.unitMomentum();
	Vector pos = part.getPosition() / Sim->units.unitLength();
	Vector oldv = pData.getOldVel() / Sim->units.unitVelocity();
	Vector newv = part.getVelocity() / Sim->units.unitVelocity();
	logfile << " delP1=" << delP.toString() << ", pos=" << pos.toString() << ", vel=" << newv.toString() << ", oldvel=" << oldv.toString() << "\n";
      }

    for (const PairEventData& pData : SDat.L2partChanges)
      {
	logfile << "    2PEvent";
      
	printData(pData.particle1_.getParticleID(),
		  pData.particle2_.getParticleID());
      
	logfile << "\n";
      }
  }

  void 
  OPTrajectory::eventUpdate(const System& sys, const NEventData& SDat, 
			    const double& dt)
  {
    logfile << std::setw(8) << Sim->eventCount
	    << " SYSTEM " << sys.getID()
	    << " TYPE " << sys.getType()
	    << " t " << Sim->systemTime / Sim->units.unitTime() 
	    << " dt " << dt / Sim->units.unitTime()
	    << "\n";

    for (const ParticleEventData& pData : SDat.L1partChanges)
      {
	const Particle& part = Sim->particles[pData.getParticleID()];
	logfile << "    1PEvent p1 " << part.getID();
	Vector delP = Sim->species[pData.getSpeciesID()]->getMass(part.getID()) * (part.getVelocity() - pData.getOldVel());
	delP /= Sim->units.unitMomentum();
	Vector pos = part.getPosition() / Sim->units.unitLength();	
	Vector oldv = pData.getOldVel() / Sim->units.unitVelocity();
	Vector newv = part.getVelocity() / Sim->units.unitVelocity();
	logfile << " delP1=" << delP.toString() << ", pos=" << pos.toString() << ", vel=" << newv.toString() << ", oldvel=" << oldv.toString() << "\n";
      }
  
    for (const PairEventData& pData : SDat.L2partChanges)
      {
	logfile << "    2PEvent";
      
	printData(pData.particle1_.getParticleID(),
		  pData.particle2_.getParticleID());
      
	logfile << "\n";
      }
  }

  void 
  OPTrajectory::output(magnet::xml::XmlStream& XML)
  {}
}
