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

const CIntEvent 
CScheduler::earliestIntEvent() const
{
#ifdef DYNAMO_DEBUG
  if (sorter->next_Data().top().type != INTERACTION)
    D_throw() << "The next event is not an Interaction event";
#endif
  
  return Sim->Dynamics.getEvent
    (Sim->vParticleList[sorter->next_ID()], 
     Sim->vParticleList[sorter->next_Data().top().p2]);
}

const CGlobEvent
CScheduler::earliestGlobEvent() const
{
#ifdef DYNAMO_DEBUG
  if (sorter->next_Data().top().type != GLOBAL)
    D_throw() << "The next event is not a Global event";
#endif

  return Sim->Dynamics.getGlobals()[sorter->next_Data().top().p2]
    ->getEvent(Sim->vParticleList[sorter->next_ID()]);
}

const CLocalEvent
CScheduler::earliestLocalEvent() const
{
#ifdef DYNAMO_DEBUG
  if (sorter->next_Data().top().type != LOCAL)
    D_throw() << "The next event is not a Local event";
#endif

  return Sim->Dynamics.getLocals()[sorter->next_Data().top().p2]
    ->getEvent(Sim->vParticleList[sorter->next_ID()]);
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

EEventType
CScheduler::nextEventType() const
{
  //Determine the next global and/or system event
  Iflt tmpt = HUGE_VAL;

  sorter->sort();

  if (!Sim->Dynamics.getSystemEvents().empty())
    tmpt =(*min_element(Sim->Dynamics.getSystemEvents().begin(),
			Sim->Dynamics.getSystemEvents().end()
			))->getdt();
  
#ifdef DYNAMO_DEBUG
  if (sorter->next_Data().empty())
    D_throw() << "Next particle list is empty but top of list!";
#endif  
  
#ifdef DYNAMO_UpdateCollDebug
  std::cerr << "\nNext eventdt = " << sorter->next_dt();
#endif



  while (sorter->next_dt() < tmpt)
    {
      //Return it if its not an INTERACTION
      if (sorter->next_Data().top().type != INTERACTION)
	return sorter->next_Data().top().type;

      //Check the INTERACTION is valid first
      if (sorter->next_Data().top().collCounter2 
	  == eventCount[sorter->next_Data().top().p2])
	return INTERACTION;

      //Not valid, update the list
#ifdef DYNAMO_UpdateCollDebug
      std::cerr << "\nEvent invalid, popping and updating" 
		<< sorter->next_dt();
#endif
      sorter->next_Data().pop();
      sorter->update(sorter->next_ID());
      sorter->sort();

#ifdef DYNAMO_DEBUG
      if (sorter->next_Data().empty())
	D_throw() << "Next particle list is empty but top of list!";
#endif  

    }
  
  //The other events didn't win
  return SYSTEM;
}
