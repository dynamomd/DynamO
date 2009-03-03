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

#include "PBCSentinel.hpp"
#include "globEvent.hpp"
#include "../NparticleEventData.hpp"
#include "../../base/is_simdata.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../../schedulers/scheduler.hpp"

CGPBCSentinel::CGPBCSentinel(DYNAMO::SimData* nSim, const std::string& name):
  CGlobal(nSim, "PBCSentinel"),
  maxintdist(0)
{
  globName = name;
  I_cout() << "PBCSentinel Loaded";
}

CGPBCSentinel::CGPBCSentinel(const XMLNode &XML, DYNAMO::SimData* ptrSim):
  CGlobal(ptrSim, "PBCSentinel"),
  maxintdist(0)
{
  operator<<(XML);

  I_cout() << "PBCSentinel Loaded";
}

void 
CGPBCSentinel::initialise(size_t nID)
{
  ID=nID;
  
  maxintdist = Sim->Dynamics.getLongestInteraction();
  
  cachedTimes.resize(Sim->lN);
  
  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    cachedTimes[part.getID()] = Sim->dSysTime;

  Sim->registerParticleUpdateFunc
    (fastdelegate::MakeDelegate(this, &CGPBCSentinel::particlesUpdated));
}

void 
CGPBCSentinel::particlesUpdated(const CNParticleData& PDat)
{
  BOOST_FOREACH(const C1ParticleData& pdat, PDat.L1partChanges)
    cachedTimes[pdat.getParticle().getID()] = Sim->dSysTime;
  
  BOOST_FOREACH(const C2ParticleData& pdat, PDat.L2partChanges)
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

CGlobEvent 
CGPBCSentinel::getEvent(const CParticle& part) const
{
  Iflt dt 
    = Sim->Dynamics.Liouvillean().getPBCSentinelTime(part, maxintdist)
    - (Sim->dSysTime - cachedTimes[part.getID()]);
 
  return CGlobEvent(part, dt, VIRTUAL, *this);
}

void 
CGPBCSentinel::runEvent(const CParticle& part) const
{
  Sim->Dynamics.Liouvillean().updateParticle(part);

  CGlobEvent iEvent(getEvent(part));

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
  
  //dynamics must be updated first
  Sim->Dynamics.stream(iEvent.getdt());

  CNParticleData EDat(C1ParticleData(part, Sim->Dynamics.getSpecies(part), VIRTUAL));
  
  //Now we're past the event, update the scheduler and plugins

  //This update actually updates the sentinel's events
  Sim->signalParticleUpdate(EDat);

  Sim->ptrScheduler->fullUpdate(part);
  
  BOOST_FOREACH(smrtPlugPtr<COutputPlugin> & Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(iEvent, EDat);
}

void 
CGPBCSentinel::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "PBCSentinel"
      << xmlw::attr("Name") << globName;
}
