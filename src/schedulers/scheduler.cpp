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

#include "scheduler.hpp"
#include "../dynamics/interactions/intEvent.hpp"
#include "../dynamics/globals/global.hpp"
#include "../dynamics/globals/globEvent.hpp"
#include "../dynamics/locals/local.hpp"
#include "../dynamics/systems/system.hpp"
#include "../base/is_simdata.hpp"
#include "../extcode/xmlwriter.hpp"
#include "../extcode/xmlParser.h"
#include "include.hpp"

CScheduler::CScheduler(const DYNAMO::SimData* const tmp, const char * aName,
		       CSSorter* nS):
  SimBase_const(tmp, aName, IC_purple),
  sorter(nS)
{}

CScheduler::~CScheduler()
{}

CScheduler* 
CScheduler::getClass(const XMLNode& XML, const DYNAMO::SimData* Sim)
{
  if (!strcmp(XML.getAttribute("Type"),"NeighbourList"))
    return new CSNeighbourList(XML,Sim);
  else if (!strcmp(XML.getAttribute("Type"),"Dumb"))
    return new CSDumb(XML,Sim);
  else 
    D_throw() << "Unknown type of Scheduler encountered";
}

xmlw::XmlStream& operator<<(xmlw::XmlStream& XML, 
			    const CScheduler& g)
{
  g.outputXML(XML);
  return XML;
}

void 
CScheduler::rebuildSystemEvents() const
{
  (*sorter)[Sim->lN].clear();

  BOOST_FOREACH(const smrtPlugPtr<CSystem>& sysptr, 
		Sim->Dynamics.getSystemEvents())
    sorter->push(intPart(sysptr->getdt(), SYSTEM, sysptr->getID(), 0), Sim->lN);

  sorter->update(Sim->lN);
}

void 
CScheduler::popNextEvent()
{
  (*sorter)[sorter->next_ID()].pop();
}

void 
CScheduler::pushEvent(const CParticle& part,
		      const intPart& newevent)
{
  sorter->push(newevent, part.getID());
}

void 
CScheduler::sort(const CParticle& part)
{
  sorter->update(part.getID());
}

void 
CScheduler::invalidateEvents(const CParticle& part)
{
  //Invalidate previous entries
  ++eventCount[part.getID()];
  (*sorter)[part.getID()].clear();
}

void
CScheduler::runNextEvent() const
{
  sorter->sort();
  
#ifdef DYNAMO_DEBUG
  if (sorter->next_Data().empty())
    D_throw() << "Next particle list is empty but top of list!";
#endif  
    
  while ((sorter->next_Data().top().type == INTERACTION)
	 && (sorter->next_Data().top().collCounter2 
	     != eventCount[sorter->next_Data().top().p2]))
    {
      //Not valid, update the list
      sorter->next_Data().pop();
      sorter->update(sorter->next_ID());
      sorter->sort();
      
#ifdef DYNAMO_DEBUG
      if (sorter->next_Data().empty())
	D_throw() << "Next particle list is empty but top of list!";
#endif  
    }
  
  switch (sorter->next_Data().top().type)
    {
    case INTERACTION:
      Sim->Dynamics.runIntEvent
	(Sim->vParticleList[sorter->next_ID()], 
	 Sim->vParticleList[sorter->next_Data().top().p2]);	  
      break;
    case GLOBAL:
      Sim->Dynamics.getGlobals()[sorter->next_Data().top().p2]
	->runEvent(Sim->vParticleList[sorter->next_ID()]);
      break;	           
    case LOCAL:
      Sim->Dynamics.getLocals()[sorter->next_Data().top().p2]
	->runEvent(Sim->vParticleList[sorter->next_ID()]);	  
      break;
    case SYSTEM:
      Sim->Dynamics.getSystemEvents()[sorter->next_Data().top().p2]
	->runEvent();
      //This saves the system events rebuilding themselves
      rebuildSystemEvents();
      break;
    default:
      D_throw() << "Unhandled event type requested to be run";
    }
}
