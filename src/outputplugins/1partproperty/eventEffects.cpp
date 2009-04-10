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

#include "eventEffects.hpp"
#include "../../base/is_simdata.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include <boost/foreach.hpp>
#include "../../dynamics/include.hpp"
#include "../../dynamics/interactions/include.hpp"
#include "../1partproperty/kenergy.hpp"
#include "../../datatypes/vector.xml.hpp"

COPEventEffects::COPEventEffects(const DYNAMO::SimData* tmp, const XMLNode&):
  COutputPlugin(tmp,"EventEffects")
{
}

void 
COPEventEffects::initialise()
{}

COPEventEffects::~COPEventEffects()
{}

void 
COPEventEffects::eventUpdate(const CIntEvent& iEvent, const C2ParticleData& Pdat)
{
  newEvent(iEvent.getType(),getClassKey(iEvent),
	   Pdat.particle1_.getDeltaKE(),
	   -Pdat.dP);

  newEvent(iEvent.getType(),getClassKey(iEvent),
	   Pdat.particle2_.getDeltaKE(),
	   Pdat.dP);
}


void 
COPEventEffects::eventUpdate(const CGlobEvent& globEvent, const CNParticleData& SDat)
{
  BOOST_FOREACH(const C1ParticleData& pData, SDat.L1partChanges)
    newEvent(globEvent.getType(),getClassKey(globEvent),
	     pData.getDeltaKE(),
	     pData.getDeltaP());
  
  BOOST_FOREACH(const C2ParticleData& pData, SDat.L2partChanges)
    {
      newEvent(globEvent.getType(),getClassKey(globEvent),
	       pData.particle1_.getDeltaKE(),
	       -pData.dP);
      
      newEvent(globEvent.getType(),getClassKey(globEvent),
	       pData.particle2_.getDeltaKE(),
	       pData.dP);
    }
}

void 
COPEventEffects::eventUpdate(const CLocalEvent& localEvent, const CNParticleData& SDat)
{
  BOOST_FOREACH(const C1ParticleData& pData, SDat.L1partChanges)
    newEvent(localEvent.getType(),getClassKey(localEvent),
	     pData.getDeltaKE(),
	     pData.getDeltaP());
  
  BOOST_FOREACH(const C2ParticleData& pData, SDat.L2partChanges)
    {
      newEvent(localEvent.getType(),getClassKey(localEvent),
	       pData.particle1_.getDeltaKE(),
	       -pData.dP);
      
      newEvent(localEvent.getType(),getClassKey(localEvent),
	       pData.particle2_.getDeltaKE(),
	       pData.dP);
    }
}

void 
COPEventEffects::eventUpdate(const CSystem& sysEvent, const CNParticleData& SDat, const Iflt&)
{
  BOOST_FOREACH(const C1ParticleData& pData, SDat.L1partChanges)
    newEvent(sysEvent.getType(),getClassKey(sysEvent),
	     pData.getDeltaKE(),
	     pData.getDeltaP());
  
  BOOST_FOREACH(const C2ParticleData& pData, SDat.L2partChanges)
    {
      newEvent(sysEvent.getType(),getClassKey(sysEvent),
	       pData.particle1_.getDeltaKE(),
	       -pData.dP);
      
      newEvent(sysEvent.getType(),getClassKey(sysEvent),
	       pData.particle2_.getDeltaKE(),
	       pData.dP);
    }
}


void 
COPEventEffects::newEvent(const EEventType& eType, const classKey& ck, 
			  const Iflt& deltaKE, const Vector & delP)
{
  counterData& ref(counters[eventKey(ck, eType)]);
  ref.energyLoss += deltaKE;
  ref.momentumChange += delP;
}

void
COPEventEffects::output(xmlw::XmlStream &XML)
{
  XML << xmlw::tag("EventEffects");

  typedef std::pair<const eventKey, counterData> locPair;  
  
  BOOST_FOREACH(const locPair& ele, counters)
    {
      XML << xmlw::tag("Count")
	  << xmlw::attr("Name") << getName(ele.first.first, Sim)
	  << xmlw::attr("Event") << CIntEvent::getCollEnumName(ele.first.second)
	  << xmlw::attr("EnergyLossRate") 
	  << ele.second.energyLoss * Sim->Dynamics.units().unitTime()
	/ (Sim->dSysTime * Sim->Dynamics.units().unitEnergy())
	  << xmlw::tag("MomentumChangeRate") 
	  << ele.second.momentumChange * Sim->Dynamics.units().unitTime()
	/ (Sim->dSysTime * Sim->Dynamics.units().unitMomentum())
	  << xmlw::endtag("MomentumChangeRate") 
	  << xmlw::endtag("Count");
    }
  
  XML << xmlw::endtag("EventEffects");
}
