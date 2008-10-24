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

#include "globalCellular.hpp"
#include "../dynamics/interactions/intEvent.hpp"
#include "../extcode/threadpool.hpp"
#include "../simulation/particle.hpp"
#include "../dynamics/dynamics.hpp"
#include "../dynamics/liouvillean/liouvillean.hpp"
#include "../dynamics/BC/BC.hpp"
#include "../dynamics/BC/LEBC.hpp"
#include "../base/is_simdata.hpp"
#include "../base/is_base.hpp"
#include "../dynamics/globals/globEvent.hpp"
#include "../dynamics/systems/system.hpp"
#include <cmath> //for huge val
#include <boost/lexical_cast.hpp>
#include "../extcode/xmlParser.h"
#include "../dynamics/globals/global.hpp"
#include "../dynamics/globals/globEvent.hpp"
#include "../dynamics/globals/gcells.hpp"
#include "../dynamics/locals/local.hpp"
#include "../dynamics/locals/localEvent.hpp"

void 
CSGlobCellular::stream(const Iflt dt)
{
  eventHeap.stream(dt);
}

const CIntEvent 
CSGlobCellular::earliestIntEvent() const
{
#ifdef DYNAMO_DEBUG
  if (eventHeap.next_Data.top().type != INTERACTION)
    D_throw() << "The next event is not an Interaction event";
#endif
  
  return Sim->Dynamics.getEvent
    (Sim->vParticleList[eventHeap.next_ID()], 
     Sim->vParticleList[eventHeap.next_Data().top().p2]);
}

void 
CSGlobCellular::operator<<(const XMLNode& XML)
{
}

const CGlobEvent
CSGlobCellular::earliestGlobEvent() const
{
#ifdef DYNAMO_DEBUG
  if (eventHeap.next_Data.top().type != GLOBAL)
    D_throw() << "The next event is not a Global event";
#endif

  return Sim->Dynamics.getGlobals()[eventHeap.next_Data().top().p2]
    ->getEvent(Sim->vParticleList[eventHeap.next_ID()]);
}

const CLocalEvent
CSGlobCellular::earliestLocalEvent() const
{
#ifdef DYNAMO_DEBUG
  if (eventHeap.next_Data.top().type != LOCAL)
    D_throw() << "The next event is not a Local event";
#endif

  return Sim->Dynamics.getLocals()[eventHeap.next_Data().top().p2]
    ->getEvent(Sim->vParticleList[eventHeap.next_ID()]);
}

void
CSGlobCellular::initialise()
{
  if (Sim->Dynamics.BCTypeTest<CRLEBC>()
      || Sim->Dynamics.BCTypeTest<CSLEBC>())
    D_throw() << "This scheduler isn't suitable for sheared systems";

  try {
    GlobCellID = Sim->Dynamics.getGlobal("Cells")->getID();
  }
  catch(std::exception& cxp)
    {
      D_throw() << "Failed while finding the cell global event.\n"
		<< "You must have a cellular grid enabled for this scheduler.\n"
		<< "You can add one using dynamod --GCells"
		<< cxp.what();
    }

  I_cout() << "Reinitialising on collision " << Sim->lNColl;

  eventHeap.clear();
  eventHeap.resize(Sim->lN);
  eventCount.clear();
  eventCount.resize(Sim->lN, 0);

  //Now initialise the interactions
  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    addNewEvents(part);
  
  eventHeap.init();

#ifndef CBT
  I_cout() << "BPQ: Number of lists " << eventHeap.NLists();
  I_cout() << "BPQ: Scale Factor " << eventHeap.scaleFactor();
#endif
}

void 
CSGlobCellular::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "GlobalCellular";
}

CSGlobCellular::CSGlobCellular(const XMLNode& XML, const DYNAMO::SimData* Sim):
  CScheduler(Sim,"GlobalCellular")
{ 
  I_cout() << "Global Cellular Algorithmn";
  operator<<(XML);
}

CSGlobCellular::CSGlobCellular(const DYNAMO::SimData* Sim):
  CScheduler(Sim,"GlobalCellular")
{ I_cout() << "Global Cellular Algorithmn"; }

void 
CSGlobCellular::rescaleTimes(Iflt scale)
{ eventHeap.rescaleTimes(scale); }


void 
CSGlobCellular::popVirtualEvent()
{
  eventHeap[eventHeap.next_ID()].pop();
}

void 
CSGlobCellular::virtualCellNewNeighbour(const CParticle& part, const CParticle& part2)
{
  CIntEvent eevent(Sim->Dynamics.getEvent(part, part2));

  if (eevent.getType() != NONE)
    eventHeap.push(intPart(eevent, eventCount[part2.getID()]), part.getID());  
}

void 
CSGlobCellular::pushAndUpdateVirtualEvent(const CParticle& part, const intPart& newevent)
{
  eventHeap.push(newevent,part.getID());
  eventHeap.update(part.getID());
}

void 
CSGlobCellular::update(const CParticle& part)
{
  //Invalidate previous entries
  ++eventCount[part.getID()];
  eventHeap[part.getID()].clear();
  addNewEvents(part);
  eventHeap.update(part.getID());
}

ENextEvent 
CSGlobCellular::nextEventType() const
{
  //Determine the next global and/or system event
  Iflt tmpt = HUGE_VAL;

  eventHeap.sort();

  if (!Sim->Dynamics.getSystemEvents().empty())
    tmpt =(*min_element(Sim->Dynamics.getSystemEvents().begin(),
			Sim->Dynamics.getSystemEvents().end()
			))->getdt();
  
#ifdef DYNAMO_DEBUG
  if (eventHeap.next_Data().empty())
    D_throw() << "Next particle list is empty but top of list!";
#endif  
  
#ifdef DYNAMO_UpdateCollDebug
  std::cerr << "\nNext eventdt = " << eventHeap.next_dt();
#endif

  while (eventHeap.next_dt() < tmpt)
    {
     switch (eventHeap.next_Data().top().type)
      {
      case INTERACTION:
	if (eventHeap.next_Data().top().collCounter2 
	    != eventCount[eventHeap.next_Data().top().p2])
	  {
#ifdef DYNAMO_UpdateCollDebug
	    std::cerr << "\nEvent invalid, popping and updating" << eventHeap.next_dt();
#endif
	    eventHeap.next_Data().pop();
	    eventHeap.update(eventHeap.next_ID());
	    break;
	  }

	return Interaction;
      case GLOBAL:
	return Global;
	break;
      case LOCAL:
	return Local;
	break;
      default:
	D_throw() << "Unknown event type!";
      }
     eventHeap.sort();

#ifdef DYNAMO_UpdateCollDebug
     std::cerr << "\nNext eventdt = " << eventHeap.next_dt();
#endif
    }

  return System;
}

void 
CSGlobCellular::addNewEvents(const CParticle& part) const
{  
  //Add the global event
  BOOST_FOREACH(const smrtPlugPtr<CGlobal>& glob, Sim->Dynamics.getGlobals())
    if (glob->isInteraction(part))
      eventHeap.push(glob->getEvent(part), part.getID());

  const CGCells& cells(*static_cast<const CGCells*>
		       (Sim->Dynamics.getGlobals()[GlobCellID].get_ptr()));

  //Add the local cell events
  BOOST_FOREACH(const size_t& id, cells.getParticleCellData(part).locals)
    {
      const smrtPlugPtr<CLocal>& local(Sim->Dynamics.getLocals()[id]);
      if (local->isInteraction(part))
	eventHeap.push(local->getEvent(part), part.getID());
    }
      

  BOOST_FOREACH(const int& nb, cells.getCellNeighbourHood(part))
    for (int next = cells.getCellLocalParticles(nb);
	 next != -1; next = cells.getParticleData(next).next)
      if (next != static_cast<int>(part.getID()))
	{
	  CIntEvent eevent(Sim->Dynamics.getEvent(part, Sim->vParticleList[next]));
	  if (eevent.getType() != NONE)
	    eventHeap.push(intPart(eevent, eventCount[next]), part.getID());
	}
}
