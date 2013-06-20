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

#include <dynamo/outputplugins/eventEffects.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/include.hpp>
#include <dynamo/interactions/include.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
  OPEventEffects::OPEventEffects(const dynamo::Simulation* tmp, const magnet::xml::Node&):
    OutputPlugin(tmp,"EventEffects")
  {}

  void 
  OPEventEffects::initialise()
  {}

  OPEventEffects::~OPEventEffects()
  {}

  void 
  OPEventEffects::eventUpdate(const IntEvent& iEvent, const PairEventData& Pdat)
  {
    const Particle& p1 = Sim->particles[Pdat.particle1_.getParticleID()];
    const Particle& p2 = Sim->particles[Pdat.particle2_.getParticleID()];

    newEvent(iEvent.getType(),getClassKey(iEvent), 0.5 * Sim->species[p1]->getMass(p1.getID()) * (p1.getVelocity().nrm2() - Pdat.particle1_.getOldVel().nrm2()), -Pdat.impulse);
    newEvent(iEvent.getType(),getClassKey(iEvent), 0.5 * Sim->species[p2]->getMass(p2.getID()) * (p2.getVelocity().nrm2() - Pdat.particle2_.getOldVel().nrm2()), Pdat.impulse);
  }

  void 
  OPEventEffects::eventUpdate(const GlobalEvent& globEvent, const NEventData& SDat)
  {
    for (const ParticleEventData& pData : SDat.L1partChanges)
      {
	const Particle& p1 = Sim->particles[pData.getParticleID()];
	const Species& sp1 = *Sim->species[pData.getSpeciesID()];
	const double m1 = sp1.getMass(p1.getID());
	const Vector dP = m1 * (p1.getVelocity() - pData.getOldVel());

	newEvent(globEvent.getType(), getClassKey(globEvent), 0.5 * m1 * (p1.getVelocity().nrm2() - pData.getOldVel().nrm2()), dP);
      }
  
    for (const PairEventData& pData : SDat.L2partChanges)
      {
	const Particle& p1 = Sim->particles[pData.particle1_.getParticleID()];
	const Particle& p2 = Sim->particles[pData.particle2_.getParticleID()];
	const double m1 = Sim->species[p1]->getMass(p1.getID());
	const double m2 = Sim->species[p2]->getMass(p2.getID());

	newEvent(globEvent.getType(), getClassKey(globEvent),
		 0.5 * m1 * (p1.getVelocity().nrm2() - pData.particle1_.getOldVel().nrm2()),
		 -pData.impulse);
      
	newEvent(globEvent.getType(), getClassKey(globEvent),
		 0.5 * m2 * (p2.getVelocity().nrm2() - pData.particle2_.getOldVel().nrm2()),
		 pData.impulse);
      }
  }

  void 
  OPEventEffects::eventUpdate(const LocalEvent& localEvent, const NEventData& SDat)
  {
    for (const ParticleEventData& pData : SDat.L1partChanges)
      {
	const Particle& p1 = Sim->particles[pData.getParticleID()];
	const double m1 = Sim->species[p1]->getMass(p1.getID());
	const Vector dP = m1 * (p1.getVelocity() - pData.getOldVel());
	
	newEvent(localEvent.getType(),getClassKey(localEvent), 0.5 * m1 * (p1.getVelocity().nrm2() - pData.getOldVel().nrm2()), dP);
      }
  
    for (const PairEventData& pData : SDat.L2partChanges)
      {
	const Particle& p1 = Sim->particles[pData.particle1_.getParticleID()];
	const Particle& p2 = Sim->particles[pData.particle2_.getParticleID()];
	const double m1 = Sim->species[p1]->getMass(p1.getID());
	const double m2 = Sim->species[p2]->getMass(p2.getID());

	newEvent(localEvent.getType(),getClassKey(localEvent),
		 0.5 * m1 * (p1.getVelocity().nrm2() - pData.particle1_.getOldVel().nrm2()),
		 -pData.impulse);
      
	newEvent(localEvent.getType(),getClassKey(localEvent),
		 0.5 * m2 * (p2.getVelocity().nrm2() - pData.particle2_.getOldVel().nrm2()),
		 pData.impulse);
      }
  }

  void 
  OPEventEffects::eventUpdate(const System& sysEvent, const NEventData& SDat, const double&)
  {
    for (const ParticleEventData& pData : SDat.L1partChanges)
      {
	const Particle& p1 = Sim->particles[pData.getParticleID()];
	const double m1 = Sim->species[p1]->getMass(p1.getID());
	const Vector dP = m1 * (p1.getVelocity() - pData.getOldVel());
	
	newEvent(sysEvent.getType(),getClassKey(sysEvent), 0.5 * m1 * (p1.getVelocity().nrm2() - pData.getOldVel().nrm2()), dP);
      }
  
    for (const PairEventData& pData : SDat.L2partChanges)
      {
	const Particle& p1 = Sim->particles[pData.particle1_.getParticleID()];
	const Particle& p2 = Sim->particles[pData.particle2_.getParticleID()];
	const double m1 = Sim->species[p1]->getMass(p1.getID());
	const double m2 = Sim->species[p2]->getMass(p2.getID());

	newEvent(sysEvent.getType(),getClassKey(sysEvent),
		 0.5 * m1 * (p1.getVelocity().nrm2() - pData.particle1_.getOldVel().nrm2()),
		 -pData.impulse);
      
	newEvent(sysEvent.getType(),getClassKey(sysEvent),
		 0.5 * m2 * (p2.getVelocity().nrm2() - pData.particle2_.getOldVel().nrm2()),
		 pData.impulse);
      }
  }


  void 
  OPEventEffects::newEvent(const EEventType& eType, const classKey& ck, 
			   const double& deltaKE, const Vector & delP)
  {
    counterData& ref(counters[eventKey(ck, eType)]);
    ref.energyLoss += deltaKE;
    ref.momentumChange += delP;
  }

  void
  OPEventEffects::output(magnet::xml::XmlStream &XML)
  {
    XML << magnet::xml::tag("EventEffects");

    typedef std::pair<const eventKey, counterData> locPair;  
  
    for (const locPair& ele : counters)
      {
	XML << magnet::xml::tag("Count")
	    << magnet::xml::attr("Name") << getName(ele.first.first, Sim)
	    << magnet::xml::attr("Event") << ele.first.second
	    << magnet::xml::attr("EnergyLossRate") 
	    << ele.second.energyLoss * Sim->units.unitTime()
	  / (Sim->systemTime * Sim->units.unitEnergy())
	    << magnet::xml::tag("MomentumChangeRate") 
	    << ele.second.momentumChange * Sim->units.unitTime()
	  / (Sim->systemTime * Sim->units.unitMomentum())
	    << magnet::xml::endtag("MomentumChangeRate") 
	    << magnet::xml::endtag("Count");
      }
  
    XML << magnet::xml::endtag("EventEffects");
  }
}
