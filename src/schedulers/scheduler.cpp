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
CScheduler::runNextEvent()
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
  
#ifdef DYNAMO_DEBUG
  if (isnan(sorter->next_Data().top().dt))
    D_throw() << "Next event time is NaN"
	      << "\nTime to event "
	      << sorter->next_Data().top().dt
	      << "\nEvent Type = " 
	      << CIntEvent::getCollEnumName(sorter->next_Data().top().type)
	      << "\nOwner Particle = " << sorter->next_ID()
	      << "\nID2 = " << sorter->next_Data().top().p2;

  if (isinf(sorter->next_Data().top().dt))
    D_throw() << "Next event time is Inf!"
	      << "\nTime to event "
	      << sorter->next_Data().top().dt
	      << "\nEvent Type = " 
	      << CIntEvent::getCollEnumName(sorter->next_Data().top().type)
	      << "\nOwner Particle = " << sorter->next_ID()
	      << "\nID2 = " << sorter->next_Data().top().p2;

  if (sorter->next_Data().top().dt + eps < 0)
    D_throw() << "Next event time is less than -" << eps
	      << "\nTime to event "
	      << sorter->next_Data().top().dt
	      << "\nEvent Type = " 
	      << CIntEvent::getCollEnumName(sorter->next_Data().top().type)
      	      << "\nOwner Particle = " << sorter->next_ID()
	      << "\nID2 = " << sorter->next_Data().top().p2;
#endif  
  
  switch (sorter->next_Data().top().type)
    {
    case INTERACTION:
      {
	const CParticle& p1(Sim->vParticleList[sorter->next_ID()]);
	const CParticle& p2(Sim->vParticleList[sorter->next_Data().top().p2]);
	
	Sim->Dynamics.getLiouvillean().updateParticlePair(p1, p2);
	
	CIntEvent Event(Sim->Dynamics.getEvent(p1, p2));
	
	if ((fabs(1.0 - fabs(Event.getdt() / sorter->next_dt())) > 0.01)
	    && (fabs(sorter->next_dt() - Event.getdt()) > 1e-10*Sim->dSysTime))
	  {
	    I_cerr() << 
	      "A recalculated event time, performed to confirm the collision time, is"
	      "\nmore than 1% different to the FEL event time. Either due to free "
	      "\nstreaming inaccuracies or . Forcing a complete recalculation for "
	      "\nthe particles involved\n\n======Event Data======\n"
		     << Event.stringData(Sim)
		     << "\n\n Systime of original event is " 
		     << (sorter->next_dt() + Sim->dSysTime) 
	      / Sim->Dynamics.units().unitTime()
		     << "\n\n Systime of recalculated event "
		     << (Event.getdt() + Sim->dSysTime) 
	      / Sim->Dynamics.units().unitTime()
		     << "\n\n Systime is "
		     << Sim->dSysTime / Sim->Dynamics.units().unitTime();
	     
	    this->fullUpdate(p1, p2);
	    return;
	  }
	
#ifdef DYNAMO_DEBUG
	if (Event.getType() == NONE)
	  D_throw() << "A NONE event has reached the top of the queue."
	    "\nThe simulation has run out of events!"
	    "\n\n======Event Data======\n"
		    << Event.stringData(Sim);	  

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
	Sim->Dynamics.stream(Event.getdt());	
	
	Event.addTime(Sim->freestreamAcc);

	Sim->freestreamAcc = 0;

	Sim->Dynamics.getInteractions()[Event.getInteractionID()]
	  ->runEvent(p1,p2,Event);

	break;
      }
    case GLOBAL:
      {
	//We don't stream the system for globals as neighbour lists
	//optimise this (they dont need it).
	
	#ifdef DYNAMO_DEBUG
	//This makes neighbour list events appear in the output plugins
	if (dynamic_cast<const CGNeighbourList*>(Sim->Dynamics.getGlobals()[sorter->next_Data().top().p2].get_ptr())
	    != NULL)
	  BOOST_FOREACH(smrtPlugPtr<COutputPlugin> & Ptr, Sim->outputPlugins)
	    Ptr->eventUpdate(CGlobEvent(Sim->vParticleList[sorter->next_ID()], 0, 
					CELL, *Sim->Dynamics.getGlobals()[sorter->next_Data().top().p2]),
			     CNParticleData(C1ParticleData(Sim->vParticleList[sorter->next_ID()],
							   Sim->Dynamics.getSpecies(Sim->vParticleList[sorter->next_ID()]),
							   CELL)));
	#endif
	
	Sim->Dynamics.getGlobals()[sorter->next_Data().top().p2]
	  ->runEvent(Sim->vParticleList[sorter->next_ID()]);       	
	break;	           
      }
    case LOCAL:
      {
	const CParticle& part(Sim->vParticleList[sorter->next_ID()]);

	Sim->Dynamics.getLiouvillean().updateParticle(part);
	CLocalEvent iEvent(Sim->Dynamics.getLocals()[sorter->next_Data().top().p2]->getEvent(part));
	
	if (iEvent.getType() == NONE)
	  D_throw() << "No global collision found\n"
		    << iEvent.stringData(Sim);
	
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
	Sim->Dynamics.stream(iEvent.getdt());
	
	iEvent.addTime(Sim->freestreamAcc);
	Sim->freestreamAcc = 0;

	Sim->Dynamics.getLocals()[sorter->next_Data().top().p2]
	->runEvent(part, iEvent);	  
	break;
      }
    case SYSTEM:
      Sim->Dynamics.getSystemEvents()[sorter->next_Data().top().p2]
	->runEvent();
      //This saves the system events rebuilding themselves
      rebuildSystemEvents();
      break;
    default:
      D_throw() << "Unhandled event type requested to be run\n"
		<< "Type is " 
		<< CIntEvent::getCollEnumName(sorter->next_Data().top().type);
    }
}

void 
CScheduler::addInteractionEvent(const CParticle& part, 
				     const size_t& id) const
{
  const CParticle& part2(Sim->vParticleList[id]);

  Sim->Dynamics.getLiouvillean().updateParticle(part2);

  const CIntEvent& eevent(Sim->Dynamics.getEvent(part, part2));

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
  if (Sim->Dynamics.getLocals()[id]->isInteraction(part))
    sorter->push(Sim->Dynamics.getLocals()[id]->getEvent(part), part.getID());  
}

