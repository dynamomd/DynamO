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

#include "scheduler.hpp"
#include "../dynamics/interactions/intEvent.hpp"
#include "../dynamics/globals/global.hpp"
#include "../dynamics/globals/globEvent.hpp"
#include "../dynamics/locals/local.hpp"
#include "../dynamics/systems/system.hpp"
#include "../dynamics/liouvillean/liouvillean.hpp"
#include "../base/is_simdata.hpp"
#include "../extcode/xmlwriter.hpp"
#include "../extcode/xmlParser.h"
#include "include.hpp"
#include "../dynamics/units/units.hpp"

#ifdef DYNAMO_DEBUG
#include "../dynamics/globals/neighbourList.hpp"
#include "../dynamics/NparticleEventData.hpp"
#endif

CScheduler::CScheduler(DYNAMO::SimData* const tmp, const char * aName,
		       CSSorter* nS):
  SimBase(tmp, aName, IC_purple),
  sorter(nS)
{}

CScheduler::~CScheduler()
{}

CScheduler* 
CScheduler::getClass(const XMLNode& XML, DYNAMO::SimData* const Sim)
{
  if (!strcmp(XML.getAttribute("Type"),"NeighbourList"))
    return new CSNeighbourList(XML, Sim);
  else if (!strcmp(XML.getAttribute("Type"),"Dumb"))
    return new CSDumb(XML, Sim);
  else if (!strcmp(XML.getAttribute("Type"),"SystemOnly"))
    return new CSSystemOnly(XML, Sim);
  else if (!strcmp(XML.getAttribute("Type"),"Complex"))
    return new CSComplex(XML, Sim);
  else 
    D_throw() << XML.getAttribute("Type")
	      << ", Unknown type of Scheduler encountered";
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
  sorter->clearPEL(Sim->lN);

  BOOST_FOREACH(const smrtPlugPtr<CSystem>& sysptr, 
		Sim->dynamics.getSystemEvents())
    sorter->push(intPart(sysptr->getdt(), SYSTEM, sysptr->getID(), 0), Sim->lN);

  sorter->update(Sim->lN);
}

void 
CScheduler::popNextEvent()
{
  sorter->popNextPELEvent(sorter->next_ID());
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
  sorter->clearPEL(part.getID());
}

void
CScheduler::runNextEvent()
{
  sorter->sort();
  
#ifdef DYNAMO_DEBUG
  if (sorter->nextPELEmpty())
    D_throw() << "Next particle list is empty but top of list!";
#endif  
    
  while ((sorter->getNextEvent().type == INTERACTION)
	 && (sorter->getNextEvent().collCounter2 
	     != eventCount[sorter->getNextEvent().p2]))
    {
      //Not valid, update the list
      sorter->popNextEvent();
      sorter->update(sorter->next_ID());
      sorter->sort();
      
#ifdef DYNAMO_DEBUG
      if (sorter->nextPELEmpty())
	D_throw() << "Next particle list is empty but top of list!";
#endif
    }
  
#ifdef DYNAMO_DEBUG
  if (isnan(sorter->getNextEvent().dt))
    D_throw() << "Next event time is NaN"
	      << "\nTime to event "
	      << sorter->getNextEvent().dt
	      << "\nEvent Type = " 
	      << sorter->getNextEvent().type
	      << "\nOwner Particle = " << sorter->next_ID()
	      << "\nID2 = " << sorter->getNextEvent().p2;

  if (isinf(sorter->getNextEvent().dt))
    D_throw() << "Next event time is Inf!"
	      << "\nTime to event "
	      << sorter->getNextEvent().dt
	      << "\nEvent Type = " 
	      << sorter->getNextEvent().type
	      << "\nOwner Particle = " << sorter->next_ID()
	      << "\nID2 = " << sorter->getNextEvent().p2;

  if (sorter->getNextEvent().dt + eps < 0)
    D_throw() << "Next event time is less than -" << eps
	      << "\nTime to event "
	      << sorter->getNextEvent().dt
	      << "\nEvent Type = " 
	      << sorter->getNextEvent().type
      	      << "\nOwner Particle = " << sorter->next_ID()
	      << "\nID2 = " << sorter->getNextEvent().p2;
#endif  
  
  switch (sorter->getNextEvent().type)
    {
    case INTERACTION:
      {
	const CParticle& p1(Sim->vParticleList[sorter->next_ID()]);
	const CParticle& p2(Sim->vParticleList[sorter->getNextEvent().p2]);

	//Copy the FEL event
	intPart FELEvent = sorter->getNextEvent();

	//Ready the next event in the FEL
	sorter->popNextEvent();
	sorter->update(sorter->next_ID());
	sorter->sort();
	
	//Now recalculate the FEL event
	Sim->dynamics.getLiouvillean().updateParticlePair(p1, p2);       
	CIntEvent Event(Sim->dynamics.getEvent(p1, p2));
	
	if (Event.getdt() > sorter->getNextEvent().dt)
	  {
	    //The next FEL event is earlier than the recalculated event
	    //Grab the next event ID's
	    const unsigned long np1 = sorter->next_ID(),
	      np2 = sorter->getNextEvent().p2;
	    
	    //Check if the next event is just another copy of this event with possibly reversed ID's
	    if ((sorter->getNextEvent().type != INTERACTION)
		|| ((p1.getID() != np1) && (p1.getID() != np2))
		|| ((p2.getID() != np1) && (p2.getID() != np2)))
	      {
		
#ifdef DYNAMO_DEBUG
		I_cerr() << "Interaction event found to occur later than the next "
		  "FEL event [" << p1.getID() << "," << p2.getID()
			 << "] (small numerical error), recalculating";
#endif		
		this->fullUpdate(p1, p2);
		return;
	      }
	    //It's just another version of this event so we can execute it
	  }

	if (Event.getType() == NONE)
	  {
#ifdef DYNAMO_DEBUG
	    I_cerr() << "Interaction event found not to occur [" << p1.getID() << "," << p2.getID()
		     << "] (possible glancing collision canceled due to numerical error)";
#endif		
	    this->fullUpdate(p1, p2);
	    return;
	  }
	
#ifdef DYNAMO_DEBUG

	if (isnan(Event.getdt()))
	  D_throw() << "A NAN Interaction collision time has been found"
		    << Event.stringData(Sim);
	
	if (Event.getdt() == HUGE_VAL)
	  D_throw() << "An infinite Interaction (not marked as NONE) collision time has been found\n"
		    << Event.stringData(Sim);
#endif
	
	//Debug section
#ifdef DYNAMO_CollDebug
	if (p2.getID() < p2.getID())
	  std::cerr << "\nsysdt " << Event.getdt() + dSysTime
		    << "  ID1 " << p1.getID() 
		    << "  ID2 " << p2.getID()
		    << "  dt " << Event.getdt()
		    << "  Type " << CIntEvent::getCollEnumName(Event.getType());
	else
	  std::cerr << "\nsysdt " << Event.getdt() + dSysTime
		    << "  ID1 " << p2().getID() 
		    << "  ID2 " << p1().getID()
		    << "  dt " << Event.getdt()
		    << "  Type " << CIntEvent::getCollEnumName(Event.getType());
#endif

	Sim->dSysTime += Event.getdt();
	
	stream(Event.getdt());
	
	//dynamics must be updated first
	Sim->dynamics.stream(Event.getdt());	
	
	Event.addTime(Sim->freestreamAcc);

	Sim->freestreamAcc = 0;

	Sim->dynamics.getInteractions()[Event.getInteractionID()]
	  ->runEvent(p1,p2,Event);

	break;
      }
    case GLOBAL:
      {
	//We don't stream the system for globals as neighbour lists
	//optimise this (they dont need it).
	Sim->dynamics.getGlobals()[sorter->getNextEvent().p2]
	  ->runEvent(Sim->vParticleList[sorter->next_ID()]);       	
	break;	           
      }
    case LOCAL:
      {
	const CParticle& part(Sim->vParticleList[sorter->next_ID()]);

	//Copy the FEL event
	intPart FELEvent = sorter->getNextEvent();

	//Ready the next event in the FEL
	sorter->popNextEvent();
	sorter->update(sorter->next_ID());
	sorter->sort();
	
	Sim->dynamics.getLiouvillean().updateParticle(part);
	CLocalEvent iEvent(Sim->dynamics.getLocals()[FELEvent.p2]->getEvent(part));

	if (iEvent.getType() == NONE)
	  {
#ifdef DYNAMO_DEBUG
	    I_cerr() << "Local event found not to occur [" << p1.getID() << "," << p2.getID()
		     << "] (possible glancing/tenuous event canceled due to numerical error)";
#endif		
	    this->fullUpdate(part);
	    return;
	  }

	if (iEvent.getdt() > sorter->getNextEvent().dt)
	  {
#ifdef DYNAMO_DEBUG 
	    I_cerr() << "Recalculated LOCAL event time is greater than the next event time, recalculating";
#endif
	    this->fullUpdate(part);
	    return;
	  }

#ifdef DYNAMO_DEBUG 
	if (isnan(iEvent.getdt()))
	  D_throw() << "A NAN Global collision time has been found\n"
		    << iEvent.stringData(Sim);
	
	if (iEvent.getdt() == HUGE_VAL)
	  D_throw() << "An infinite (not marked as NONE) Global collision time has been found\n"
		    << iEvent.stringData(Sim);
#endif
	
	Sim->dSysTime += iEvent.getdt();
	
	stream(iEvent.getdt());
	
	//dynamics must be updated first
	Sim->dynamics.stream(iEvent.getdt());
	
	iEvent.addTime(Sim->freestreamAcc);
	Sim->freestreamAcc = 0;

	Sim->dynamics.getLocals()[FELEvent.p2]->runEvent(part, iEvent);	  
	break;
      }
    case SYSTEM:
      {
	Sim->dynamics.getSystemEvents()[sorter->getNextEvent().p2]
	  ->runEvent();
	//This saves the system events rebuilding themselves
	rebuildSystemEvents();
	break;
      }
    case VIRTUAL:
      {
	//Just recalc the events for these particles, no free
	//streaming (PBCSentinel will free stream virtual events but
	//for a specific reason)
	//I_cerr() << "VIRTUAL for " << sorter->next_ID();

	this->fullUpdate(Sim->vParticleList[sorter->next_ID()]);
	break;
      }
    case NONE:
      {
	D_throw() << "A NONE event has reached the top of the queue."
	  "\nThe simulation has run out of events! Aborting!";
      }
    default:
      D_throw() << "Unhandled event type requested to be run\n"
		<< "Type is " 
		<< sorter->getNextEvent().type
		<< "";
    }
}

void 
CScheduler::addInteractionEvent(const CParticle& part, 
				     const size_t& id) const
{
  const CParticle& part2(Sim->vParticleList[id]);

  Sim->dynamics.getLiouvillean().updateParticle(part2);

  const CIntEvent& eevent(Sim->dynamics.getEvent(part, part2));

  if (eevent.getType() != NONE)
    sorter->push(intPart(eevent, eventCount[id]), part.getID());
}

void 
CScheduler::addInteractionEventInit(const CParticle& part, 
					 const size_t& id) const
{
  //We'll be smart about memory and add evenly on initialisation. Not
  //using sorting only as it's unbalanced on a system where the
  //positions and ID's are correlated, e.g a lattice thats frozen on
  //initialisation.

  if (part.getID() % 2)
    {
      if (id % 2)
	//1st odd  2nd odd
	//Only take half these matches
	if (part.getID() > id) return;
      else
	//1st odd  2nd even
	//We allow these
	{}
    }
  else
    {
      if (id % 2)
	//1st even 2nd odd
	//As we allow odd,even we deny even,odd
	return;
      else
	//1st even 2nd even
	//No reason to use < or > but we switch it from odd,odd anyway
	if (part.getID() < id) return;
    }

  addInteractionEvent(part, id);
}

void 
CScheduler::addLocalEvent(const CParticle& part, 
			  const size_t& id) const
{
  if (Sim->dynamics.getLocals()[id]->isInteraction(part))
    sorter->push(Sim->dynamics.getLocals()[id]->getEvent(part), part.getID());  
}

