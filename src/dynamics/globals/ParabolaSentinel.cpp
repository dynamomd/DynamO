/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "ParabolaSentinel.hpp"
#include "globEvent.hpp"
#include "../NparticleEventData.hpp"
#include "../../base/is_simdata.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../../schedulers/scheduler.hpp"

CGParabolaSentinel::CGParabolaSentinel(DYNAMO::SimData* nSim, const std::string& name):
  CGlobal(nSim, "ParabolaSentinel")
{
  globName = name;
  I_cout() << "ParabolaSentinel Loaded";
}

CGParabolaSentinel::CGParabolaSentinel(const XMLNode &XML, DYNAMO::SimData* ptrSim):
  CGlobal(ptrSim, "ParabolaSentinel")
{
  operator<<(XML);

  I_cout() << "ParabolaSentinel Loaded";
}

void 
CGParabolaSentinel::initialise(size_t nID)
{
  ID=nID;
}

void 
CGParabolaSentinel::operator<<(const XMLNode& XML)
{
  try {
    globName = XML.getAttribute("Name");	
  }
  catch(...)
    {
      D_throw() << "Error loading CGParabolaSentinel";
    }
}

CGlobEvent 
CGParabolaSentinel::getEvent(const CParticle& part) const
{
  return CGlobEvent(part, Sim->dynamics.getLiouvillean().getParabolaSentinelTime(part), VIRTUAL, *this);
}

void 
CGParabolaSentinel::runEvent(const CParticle& part) const
{
  Sim->dynamics.getLiouvillean().updateParticle(part);

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
  
  Sim->dynamics.stream(iEvent.getdt());
  
  Sim->freestreamAcc += iEvent.getdt();

  Sim->ptrScheduler->fullUpdate(part);
}

void 
CGParabolaSentinel::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "ParabolaSentinel"
      << xmlw::attr("Name") << globName;
}
