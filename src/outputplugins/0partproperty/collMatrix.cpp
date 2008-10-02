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
  p2time.resize(Sim->lN, p2timeData(0.0, std::pair<unsigned int, EEventType>(0,NONE)));
}

void 
COPCollMatrix::initialise()
{
  Sim->getOutputPlugin<COPKEnergy>();
}

COPCollMatrix::~COPCollMatrix()
{}

size_t 
COPCollMatrix::getID(const CInteraction& E) const
{
  std::map<const std::string, unsigned int>::const_iterator 
    iPtr = intLookup.find(E.getName());

  if (iPtr == intLookup.end())
    {
      intLookup.insert(std::make_pair(E.getName(), ++IDcounter));
      return IDcounter;
    }
  else
    return iPtr->second;
}

size_t 
COPCollMatrix::getID(const CGlobal& E) const
{
  std::map<const std::string, unsigned int>::const_iterator 
    iPtr = globLookup.find(E.getName());
  
  if (iPtr == globLookup.end())
    {
      globLookup.insert(std::make_pair(E.getName(), ++IDcounter));
      return IDcounter;
    }
  else
    return iPtr->second;
}

size_t 
COPCollMatrix::getID(const CSystem& E) const
{
  std::map<const std::string, unsigned int>::const_iterator 
    iPtr = sysLookup.find(E.getName());
  
  if (iPtr == sysLookup.end())
    {
      sysLookup.insert(std::make_pair(E.getName(), ++IDcounter));
      return IDcounter;
    }
  else
    return iPtr->second;
}

void 
COPCollMatrix::eventUpdate(const CIntEvent& iEvent, const C2ParticleData&)
{
  newEvent(getID(iEvent.getInteraction()), iEvent.getParticle1(), iEvent.getType());
  newEvent(getID(iEvent.getInteraction()), iEvent.getParticle2(), iEvent.getType());
}

void 
COPCollMatrix::eventUpdate(const CGlobEvent& globEvent, const CNParticleData& SDat)
{
  BOOST_FOREACH(const C1ParticleData& pData, SDat.L1partChanges)
    newEvent(getID(globEvent.getGlobal()),pData.getParticle(), pData.getType());  
  
  BOOST_FOREACH(const C2ParticleData& pData, SDat.L2partChanges)
    {
      newEvent(getID(globEvent.getGlobal()), pData.particle1_.getParticle(), pData.getType());  
      newEvent(getID(globEvent.getGlobal()), pData.particle2_.getParticle(), pData.getType());  
    }
}

void 
COPCollMatrix::eventUpdate(const CSystem& sysEvent, const CNParticleData& SDat, const Iflt&)
{
  BOOST_FOREACH(const C1ParticleData& pData, SDat.L1partChanges)
    newEvent(getID(sysEvent), pData.getParticle(), pData.getType());  
  
  BOOST_FOREACH(const C2ParticleData& pData, SDat.L2partChanges)
    {
      newEvent(getID(sysEvent), pData.particle1_.getParticle(), pData.getType());  
      newEvent(getID(sysEvent), pData.particle2_.getParticle(), pData.getType());  
    } 
}

void 
COPCollMatrix::newEvent(const unsigned int ID, const CParticle& part, EEventType etype)
{
  p2timeData& reftData = p2time[part.getID()];
  if (p2time[part.getID()].second.first != 0)
    {
      counterData& refCount = counters[counterKey(std::pair<unsigned int, EEventType> (ID, etype),
						  std::pair<unsigned int, EEventType> (reftData.second.first, 
										       reftData.second.second))];
      
      refCount.totalTime += Sim->dSysTime - reftData.first;
      refCount.count++;
      totalCount++;
    }
  
  reftData.first = Sim->dSysTime;
  reftData.second = std::pair<unsigned int, EEventType>(ID, etype);
}

std::string 
COPCollMatrix::getName(const unsigned int& ID) const
{
  typedef std::pair<const std::string, unsigned int> localpair;

  BOOST_FOREACH(localpair& p,intLookup)
    if (p.second == ID)
      return p.first;

  BOOST_FOREACH(localpair& p, globLookup)
    if (p.second == ID)
      return p.first;

  BOOST_FOREACH(localpair& p, sysLookup)
    if (p.second == ID)
      return p.first;

  I_throw() << "Cannot find the name for ID " << ID;
}


void
COPCollMatrix::output(xmlw::XmlStream &XML)
{
  
  XML << xmlw::tag("CollCounters") 
      << xmlw::tag("CollMatrix");  

  std::map<std::pair<unsigned int, EEventType>, unsigned long> totmap;
  
  BOOST_FOREACH(const counterElement& ele, counters)
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
      totmap[std::pair<unsigned int, EEventType>
	     (ele.first.first.first,ele.first.first.second)] += ele.second.count;
    }

  XML << xmlw::endtag("CollMatrix")
      << xmlw::tag("Totals");
  
  typedef std::pair<std::pair<const unsigned int, EEventType>, unsigned long> mappair;
  BOOST_FOREACH(const mappair& mp1, totmap)
    XML << xmlw::tag("TotCount")
	<< xmlw::attr("Name") << getName(mp1.first.first)
	<< xmlw::attr("Event") << CIntEvent::getCollEnumName(mp1.first.second)
	<< xmlw::attr("Percent") << 100.0 * ((Iflt) mp1.second) / ((Iflt) totalCount)
	<< xmlw::attr("Count") << mp1.second
	<< xmlw::attr("MFT") << (Sim->dSysTime * Sim->lN) / mp1.second
	<< xmlw::endtag("TotCount");

  XML << xmlw::endtag("Totals")
      << xmlw::endtag("CollCounters");

  double pressure = 0.0;
  BOOST_FOREACH(const mappair& mp1, totmap)
    if (mp1.first.second == CORE)
      //Check its an interaction
      try {
	const smrtPlugPtr<CInteraction>& intPtr(Sim->Dynamics.getInteraction(getName(mp1.first.first)));
	//This currently works for all interactions
	pressure += intPtr->hardCoreDiam() * mp1.second;
      }
      catch (DYNAMO::Exception &)
	{}
    else if (mp1.first.second == BOUNCE)
      try {
	const smrtPlugPtr<CInteraction>& intPtr(Sim->Dynamics.getInteraction(getName(mp1.first.first)));
	if ((dynamic_cast<const CISquareBond*>(intPtr.get_ptr()) != NULL)
	    || (dynamic_cast<const CISquareWell*>(intPtr.get_ptr()) != NULL)
	    || (dynamic_cast<const CISWSequence*>(intPtr.get_ptr()) != NULL))
	  pressure -= intPtr->maxIntDist() * mp1.second;
	else
	  I_throw() << "Could turn the BOUNCE interaction into a pressure!";
      }
      catch (DYNAMO::Exception&)
	{}

  pressure *= sqrt(PI / Sim->getOutputPlugin<COPKEnergy>()->getAvgkT()) / (2.0 * NDIM * Sim->lN * Sim->dSysTime);
  pressure += 1.0;
  //Now we have the compressibility, turn it into a pressure
  pressure *= Sim->lN * Sim->getOutputPlugin<COPKEnergy>()->getAvgkT() / Sim->Dynamics.units().simVolume();

  XML << xmlw::tag("IntPressure") << xmlw::attr("val")
      << pressure / Sim->Dynamics.units().unitPressure()
      << xmlw::endtag("IntPressure");
  

}
