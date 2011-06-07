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

#include "PBCSentinel.hpp"
#include "globEvent.hpp"
#include "../NparticleEventData.hpp"
#include "../../base/is_simdata.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../../schedulers/scheduler.hpp"
#include <magnet/xmlreader.hpp>

CGPBCSentinel::CGPBCSentinel(dynamo::SimData* nSim, const std::string& name):
  Global(nSim, "PBCSentinel"),
  maxintdist(0)
{
  globName = name;
  dout << "PBCSentinel Loaded" << std::endl;
}

CGPBCSentinel::CGPBCSentinel(const magnet::xml::Node& XML, dynamo::SimData* ptrSim):
  Global(ptrSim, "PBCSentinel"),
  maxintdist(0)
{
  operator<<(XML);

  dout << "PBCSentinel Loaded" << std::endl;
}

void 
CGPBCSentinel::initialise(size_t nID)
{
  ID=nID;
  
  maxintdist = Sim->dynamics.getLongestInteraction();
  
  cachedTimes.resize(Sim->N);
  
  BOOST_FOREACH(const Particle& part, Sim->particleList)
    {
      Sim->dynamics.getLiouvillean().updateParticle(part);

      cachedTimes[part.getID()]
	= Sim->dSysTime + Sim->dynamics.getLiouvillean().getPBCSentinelTime(part, maxintdist);
    }

  Sim->registerParticleUpdateFunc
    (magnet::function::MakeDelegate(this, &CGPBCSentinel::particlesUpdated));

}

void 
CGPBCSentinel::particlesUpdated(const NEventData& PDat)
{
  BOOST_FOREACH(const ParticleEventData& pdat, PDat.L1partChanges)
    {
      cachedTimes[pdat.getParticle().getID()] 
	= Sim->dSysTime + Sim->dynamics.getLiouvillean().getPBCSentinelTime(pdat.getParticle(), maxintdist);
    }
  
  BOOST_FOREACH(const PairEventData& pdat, PDat.L2partChanges)
    {
      cachedTimes[pdat.particle1_.getParticle().getID()]
	= Sim->dSysTime 
	+ Sim->dynamics.getLiouvillean().getPBCSentinelTime(pdat.particle1_.getParticle(), maxintdist);

      cachedTimes[pdat.particle2_.getParticle().getID()]
	= Sim->dSysTime
	+ Sim->dynamics.getLiouvillean().getPBCSentinelTime(pdat.particle2_.getParticle(), maxintdist);
    }
}

void 
CGPBCSentinel::operator<<(const magnet::xml::Node& XML)
{
  globName = XML.getAttribute("Name");	
}

GlobalEvent 
CGPBCSentinel::getEvent(const Particle& part) const
{
  return GlobalEvent(part, cachedTimes[part.getID()] - Sim->dSysTime, VIRTUAL, *this);
}

void 
CGPBCSentinel::runEvent(const Particle& part, const double) const
{
  GlobalEvent iEvent(getEvent(part));

#ifdef DYNAMO_DEBUG 
  if (boost::math::isnan(iEvent.getdt()))
    M_throw() << "A NAN Interaction collision time has been found"
	      << iEvent.stringData(Sim);
  
  if (iEvent.getdt() == HUGE_VAL)
    M_throw() << "An infinite Interaction (not marked as NONE) collision time has been found\n"
	      << iEvent.stringData(Sim);
#endif

  Sim->dSysTime += iEvent.getdt();
    
  Sim->ptrScheduler->stream(iEvent.getdt());
  
  Sim->dynamics.stream(iEvent.getdt());

  Sim->dynamics.getLiouvillean().updateParticle(part);
  cachedTimes[part.getID()] 
    = Sim->dSysTime + Sim->dynamics.getLiouvillean().getPBCSentinelTime(part, maxintdist);

#ifdef DYNAMO_DEBUG
  iEvent.addTime(Sim->freestreamAcc);
  
  Sim->freestreamAcc = 0;

  NEventData EDat(ParticleEventData(part, Sim->dynamics.getSpecies(part), VIRTUAL));

  Sim->signalParticleUpdate(EDat);

  BOOST_FOREACH(magnet::ClonePtr<OutputPlugin> & Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(iEvent, EDat);
#else
  Sim->freestreamAcc += iEvent.getdt();
#endif

  Sim->ptrScheduler->fullUpdate(part);
}

void 
CGPBCSentinel::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Type") << "PBCSentinel"
      << xml::attr("Name") << globName;
}
