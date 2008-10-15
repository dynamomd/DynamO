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

#include "collMatrix.hpp"
#include "../../base/is_simdata.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include <boost/foreach.hpp>
#include "../../dynamics/include.hpp"
#include "../../dynamics/interactions/include.hpp"
#include "../1partproperty/kenergy.hpp"


COPCollMatrix::COPCollMatrix(const DYNAMO::SimData* tmp):
  COutputPlugin(tmp,"CollisionMatrix"),
  totalCount(0),
  IDcounter(0)
{
  lastEvent.resize(Sim->lN, lastEventData(Sim->dSysTime, eventKey(classKey(0, NOEventClass), NONE)));
}

void 
COPCollMatrix::initialise()
{}

COPCollMatrix::~COPCollMatrix()
{}

void 
COPCollMatrix::eventUpdate(const CIntEvent& iEvent, const C2ParticleData&)
{
  newEvent(iEvent.getParticle1(), iEvent.getType(), 
	   iEvent.getInteraction().getID(), InteractionClass);

  newEvent(iEvent.getParticle2(), iEvent.getType(), 
	   iEvent.getInteraction().getID(), InteractionClass);
}


void 
COPCollMatrix::eventUpdate(const CGlobEvent& globEvent, const CNParticleData& SDat)
{
  BOOST_FOREACH(const C1ParticleData& pData, SDat.L1partChanges)
    newEvent(pData.getParticle(), pData.getType(), globEvent.getGlobal().getID(), GlobalClass);  
  
  BOOST_FOREACH(const C2ParticleData& pData, SDat.L2partChanges)
    {
      newEvent(pData.particle1_.getParticle(), pData.getType(), globEvent.getGlobal().getID(), GlobalClass);  
      newEvent(pData.particle2_.getParticle(), pData.getType(), globEvent.getGlobal().getID(), GlobalClass);  
    }
}

void 
COPCollMatrix::eventUpdate(const CSystem& sysEvent, const CNParticleData& SDat, const Iflt&)
{
  BOOST_FOREACH(const C1ParticleData& pData, SDat.L1partChanges)
    newEvent(pData.getParticle(), pData.getType(), sysEvent.getID(), SystemClass);  
  
  BOOST_FOREACH(const C2ParticleData& pData, SDat.L2partChanges)
    {
      newEvent(pData.particle1_.getParticle(), pData.getType(), sysEvent.getID(), SystemClass);  
      newEvent(pData.particle2_.getParticle(), pData.getType(), sysEvent.getID(), SystemClass);  
    } 
}


void 
COPCollMatrix::newEvent(const CParticle& part, const EEventType& etype, const size_t& ID, const eventClass& eclass)
{
  if (lastEvent[part.getID()].second.first.second != NOEventClass)
    {      
      counterData& refCount = counters[counterKey(eventKey(classKey(ID, eclass), etype),
						  lastEvent[part.getID()].second)];
      
      refCount.totalTime += Sim->dSysTime - lastEvent[part.getID()].first;
      ++(refCount.count);
      ++(totalCount);
    }
  
  lastEvent[part.getID()].first = Sim->dSysTime;
  lastEvent[part.getID()].second = eventKey(classKey(ID, eclass), etype);
}

std::string 
COPCollMatrix::getName(const classKey& key) const
{
  switch (key.second)
    {
    case InteractionClass:
      return Sim->Dynamics.getInteractions()[key.first]->getName();
      break;
    case GlobalClass:
      return Sim->Dynamics.getGlobals()[key.first]->getName();
      break;
    case SystemClass:
      return Sim->Dynamics.getSystemEvents()[key.first]->getName();
      break;
    default:
      D_throw() << "Collision matrix found an unknown event class";
    }

}

void
COPCollMatrix::output(xmlw::XmlStream &XML)
{
  
  XML << xmlw::tag("CollCounters") 
      << xmlw::tag("CollMatrix");  
  
  std::map<eventKey, std::pair<unsigned long long, Iflt> > totmap;
  
  typedef std::pair<const counterKey, counterData> locPair;

  BOOST_FOREACH(const locPair& ele, counters)
    {
      XML << xmlw::tag("Count")
	  << xmlw::attr("Event") << CIntEvent::getCollEnumName(ele.first.first.second)
	  << xmlw::attr("Name") << getName(ele.first.first.first)
	  << xmlw::attr("lastEvent") << CIntEvent::getCollEnumName(ele.first.second.second)
	  << xmlw::attr("lastName") << getName(ele.first.second.first)
	  << xmlw::attr("Percent") << 100.0 * ((Iflt) ele.second.count) 
	/ ((Iflt) totalCount)
	  << xmlw::attr("mft") << ele.second.totalTime
	/ (Sim->Dynamics.units().unitTime() * ((Iflt) ele.second.count))
	  << xmlw::endtag("Count");
      
      //Add the total count
      totmap[ele.first.first].first += ele.second.count;

      //Add the rate
      totmap[ele.first.first].second += ((Iflt) ele.second.count) 
	/ ele.second.totalTime;
    }

  XML << xmlw::endtag("CollMatrix")
      << xmlw::tag("Totals");
  
  typedef std::pair<eventKey, std::pair<unsigned long long, Iflt> > mappair;

  BOOST_FOREACH(const mappair& mp1, totmap)
    XML << xmlw::tag("TotCount")
	<< xmlw::attr("Name") << getName(mp1.first.first)
	<< xmlw::attr("Event") << CIntEvent::getCollEnumName(mp1.first.second)
	<< xmlw::attr("Percent") << 100.0 * ((Iflt) mp1.second.first) / ((Iflt) totalCount)
	<< xmlw::attr("Count") << mp1.second.first
	<< xmlw::attr("MFT") << 1.0 / 
    (mp1.second.second * Sim->Dynamics.units().unitTime())
	<< xmlw::endtag("TotCount");

  XML << xmlw::endtag("Totals")
      << xmlw::endtag("CollCounters");
}
