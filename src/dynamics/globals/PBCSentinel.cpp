/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

CGPBCSentinel::CGPBCSentinel(DYNAMO::SimData* nSim, const std::string& name):
  Global(nSim, "PBCSentinel"),
  maxintdist(0)
{
  globName = name;
  I_cout() << "PBCSentinel Loaded";
}

CGPBCSentinel::CGPBCSentinel(const XMLNode &XML, DYNAMO::SimData* ptrSim):
  Global(ptrSim, "PBCSentinel"),
  maxintdist(0)
{
  operator<<(XML);

  I_cout() << "PBCSentinel Loaded";
}

void 
CGPBCSentinel::initialise(size_t nID)
{
  ID=nID;
  
  maxintdist = Sim->dynamics.getLongestInteraction();
  
  cachedTimes.resize(Sim->N);
  
  BOOST_FOREACH(const Particle& part, Sim->particleList)
    cachedTimes[part.getID()] = Sim->dSysTime;

  Sim->registerParticleUpdateFunc
    (fastdelegate::MakeDelegate(this, &CGPBCSentinel::particlesUpdated));

}

void 
CGPBCSentinel::particlesUpdated(const NEventData& PDat)
{
  BOOST_FOREACH(const ParticleEventData& pdat, PDat.L1partChanges)
    cachedTimes[pdat.getParticle().getID()] = Sim->dSysTime;
  
  BOOST_FOREACH(const PairEventData& pdat, PDat.L2partChanges)
    {
      cachedTimes[pdat.particle1_.getParticle().getID()] = Sim->dSysTime;
      cachedTimes[pdat.particle2_.getParticle().getID()] = Sim->dSysTime;
    }
}

void 
CGPBCSentinel::operator<<(const XMLNode& XML)
{
  try {
    globName = XML.getAttribute("Name");	
  }
  catch(...)
    {
      D_throw() << "Error loading CGPBCSentinel";
    }
}

GlobalEvent 
CGPBCSentinel::getEvent(const Particle& part) const
{
  Iflt dt 
    = Sim->dynamics.getLiouvillean().getPBCSentinelTime(part, maxintdist)
    - (Sim->dSysTime - cachedTimes[part.getID()]);
 
  return GlobalEvent(part, dt, VIRTUAL, *this);
}

void 
CGPBCSentinel::runEvent(const Particle& part) const
{
  Sim->dynamics.getLiouvillean().updateParticle(part);

  GlobalEvent iEvent(getEvent(part));

#ifdef DYNAMO_DEBUG 
  if (isnan(iEvent.getdt()))
    D_throw() << "A NAN Interaction collision time has been found"
	      << iEvent.stringData(Sim);
  
  if (iEvent.getdt() == HUGE_VAL)
    D_throw() << "An infinite Interaction (not marked as NONE) collision time has been found\n"
	      << iEvent.stringData(Sim);
#endif

  Sim->dSysTime += iEvent.getdt();
    
  Sim->ptrScheduler->stream(iEvent.getdt());
  
  Sim->dynamics.stream(iEvent.getdt());

  cachedTimes[part.getID()] = Sim->dSysTime;

  Sim->freestreamAcc += iEvent.getdt();

  Sim->ptrScheduler->fullUpdate(part);
}

void 
CGPBCSentinel::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "PBCSentinel"
      << xmlw::attr("Name") << globName;
}
