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

#include <dynamo/outputplugins/1partproperty/eventEffects.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/include.hpp>
#include <dynamo/interactions/include.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>

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
    newEvent(iEvent.getType(),getClassKey(iEvent),
	     Pdat.particle1_.getDeltaKE(),
	     -Pdat.impulse);

    newEvent(iEvent.getType(),getClassKey(iEvent),
	     Pdat.particle2_.getDeltaKE(),
	     Pdat.impulse);
  }

  void 
  OPEventEffects::eventUpdate(const GlobalEvent& globEvent, const NEventData& SDat)
  {
    BOOST_FOREACH(const ParticleEventData& pData, SDat.L1partChanges)
      {
	const Particle& p1 = Sim->particles[pData.getParticleID()];
	const Species& sp1 = *Sim->species[pData.getSpeciesID()];
	
	Vector dP = sp1.getMass(p1.getID()) * (p1.getVelocity() - pData.getOldVel());

	newEvent(globEvent.getType(),getClassKey(globEvent),
		 pData.getDeltaKE(), dP);
      }
  
    BOOST_FOREACH(const PairEventData& pData, SDat.L2partChanges)
      {
	newEvent(globEvent.getType(),getClassKey(globEvent),
		 pData.particle1_.getDeltaKE(),
		 -pData.impulse);
      
	newEvent(globEvent.getType(),getClassKey(globEvent),
		 pData.particle2_.getDeltaKE(),
		 pData.impulse);
      }
  }

  void 
  OPEventEffects::eventUpdate(const LocalEvent& localEvent, const NEventData& SDat)
  {
    BOOST_FOREACH(const ParticleEventData& pData, SDat.L1partChanges)
      {
	const Particle& p1 = Sim->particles[pData.getParticleID()];
	const Species& sp1 = *Sim->species[pData.getSpeciesID()];
	
	Vector dP = sp1.getMass(p1.getID()) * (p1.getVelocity() - pData.getOldVel());
	
	newEvent(localEvent.getType(),getClassKey(localEvent),
		 pData.getDeltaKE(), dP);
      }
  
    BOOST_FOREACH(const PairEventData& pData, SDat.L2partChanges)
      {
	newEvent(localEvent.getType(),getClassKey(localEvent),
		 pData.particle1_.getDeltaKE(),
		 -pData.impulse);
      
	newEvent(localEvent.getType(),getClassKey(localEvent),
		 pData.particle2_.getDeltaKE(),
		 pData.impulse);
      }
  }

  void 
  OPEventEffects::eventUpdate(const System& sysEvent, const NEventData& SDat, const double&)
  {
    BOOST_FOREACH(const ParticleEventData& pData, SDat.L1partChanges)
      {
	const Particle& p1 = Sim->particles[pData.getParticleID()];
	const Species& sp1 = *Sim->species[pData.getSpeciesID()];
	
	Vector dP = sp1.getMass(p1.getID()) * (p1.getVelocity() - pData.getOldVel());
	
	newEvent(sysEvent.getType(),getClassKey(sysEvent),
		 pData.getDeltaKE(), dP);
      }
  
    BOOST_FOREACH(const PairEventData& pData, SDat.L2partChanges)
      {
	newEvent(sysEvent.getType(),getClassKey(sysEvent),
		 pData.particle1_.getDeltaKE(),
		 -pData.impulse);
      
	newEvent(sysEvent.getType(),getClassKey(sysEvent),
		 pData.particle2_.getDeltaKE(),
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
  
    BOOST_FOREACH(const locPair& ele, counters)
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
