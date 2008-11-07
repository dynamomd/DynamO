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

#include "neighbourlist.hpp"
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
#include "../dynamics/globals/neighbourList.hpp"
#include "../dynamics/locals/local.hpp"
#include "../dynamics/locals/localEvent.hpp"
#include <functional>

void 
CSNeighbourList::stream(const Iflt& dt)
{
  sorter->stream(dt);
}

const CIntEvent 
CSNeighbourList::earliestIntEvent() const
{
#ifdef DYNAMO_DEBUG
  if (sorter->next_Data().top().type != INTERACTION)
    D_throw() << "The next event is not an Interaction event";
#endif
  
  return Sim->Dynamics.getEvent
    (Sim->vParticleList[sorter->next_ID()], 
     Sim->vParticleList[sorter->next_Data().top().p2]);
}

void 
CSNeighbourList::operator<<(const XMLNode& XML)
{
}

const CGlobEvent
CSNeighbourList::earliestGlobEvent() const
{
#ifdef DYNAMO_DEBUG
  if (sorter->next_Data().top().type != GLOBAL)
    D_throw() << "The next event is not a Global event";
#endif

  return Sim->Dynamics.getGlobals()[sorter->next_Data().top().p2]
    ->getEvent(Sim->vParticleList[sorter->next_ID()]);
}

const CLocalEvent
CSNeighbourList::earliestLocalEvent() const
{
#ifdef DYNAMO_DEBUG
  if (sorter->next_Data().top().type != LOCAL)
    D_throw() << "The next event is not a Local event";
#endif

  return Sim->Dynamics.getLocals()[sorter->next_Data().top().p2]
    ->getEvent(Sim->vParticleList[sorter->next_ID()]);
}

void
CSNeighbourList::initialise()
{
  if (Sim->Dynamics.BCTypeTest<CRLEBC>()
      || Sim->Dynamics.BCTypeTest<CSLEBC>())
    D_throw() << "This scheduler isn't suitable for sheared systems";
  
  try {
    NBListID = Sim->Dynamics.getGlobal("SchedulerNBList")->getID();
  }
  catch(std::exception& cxp)
    {
      D_throw() << "Failed while finding the neighbour list global.\n"
		<< "You must have a neighbour list enabled for this\n"
		<< "scheduler called SchedulerNBList.\n"
		<< cxp.what();
    }

  I_cout() << "Reinitialising on collision " << Sim->lNColl;

  sorter->clear();
  sorter->resize(Sim->lN);
  eventCount.clear();
  eventCount.resize(Sim->lN, 0);

  //Now initialise the interactions
  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    addNewEvents(part);
  
  sorter->init();
  
  //Register the new neighbour function with the cellular tracker
  if (!cellChange.connected())
    cellChange 
      = static_cast<const CGNeighbourList&>(*(Sim->Dynamics.getGlobals()[NBListID]))
      .registerCellTransitionNewNeighbourCallBack
      (boost::bind(&CSNeighbourList::virtualCellNewNeighbour, this, _1, _2));
  
  if (!reinit.connected())
    reinit 
      = static_cast<const CGNeighbourList&>(*(Sim->Dynamics.getGlobals()[NBListID]))
      .registerReInitNotify
      (boost::bind(&CSNeighbourList::initialise, this));
}

void 
CSNeighbourList::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "GlobalCellular";
}

CSNeighbourList::CSNeighbourList(const XMLNode& XML, const DYNAMO::SimData* Sim):
  CScheduler(Sim,"GlobalCellular")
{ 
  I_cout() << "Global Cellular Algorithmn";
  operator<<(XML);
}

CSNeighbourList::CSNeighbourList(const DYNAMO::SimData* Sim):
  CScheduler(Sim,"GlobalCellular")
{ I_cout() << "Global Cellular Algorithmn"; }

void 
CSNeighbourList::rescaleTimes(const Iflt& scale)
{ sorter->rescaleTimes(scale); }


void 
CSNeighbourList::popVirtualEvent()
{
  (*sorter)[sorter->next_ID()].pop();
}

void 
CSNeighbourList::virtualCellNewNeighbour(const CParticle& part, const CParticle& part2)
{
  CIntEvent eevent(Sim->Dynamics.getEvent(part, part2));

  if (eevent.getType() != NONE)
    sorter->push(intPart(eevent, eventCount[part2.getID()]), part.getID());  
}

void 
CSNeighbourList::pushAndUpdateVirtualEvent(const CParticle& part, const intPart& newevent)
{
  sorter->push(newevent,part.getID());
  sorter->update(part.getID());
}

void 
CSNeighbourList::update(const CParticle& part)
{
  //Invalidate previous entries
  ++eventCount[part.getID()];
  (*sorter)[part.getID()].clear();
  addNewEvents(part);
  sorter->update(part.getID());
}

ENextEvent 
CSNeighbourList::nextEventType() const
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
     switch (sorter->next_Data().top().type)
      {
      case INTERACTION:
	if (sorter->next_Data().top().collCounter2 
	    != eventCount[sorter->next_Data().top().p2])
	  {
#ifdef DYNAMO_UpdateCollDebug
	    std::cerr << "\nEvent invalid, popping and updating" << sorter->next_dt();
#endif
	    sorter->next_Data().pop();
	    sorter->update(sorter->next_ID());
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
     sorter->sort();

#ifdef DYNAMO_UpdateCollDebug
     std::cerr << "\nNext eventdt = " << sorter->next_dt();
#endif
    }

  return System;
}

void 
CSNeighbourList::addInteractionEvent(const CParticle& part, 
				    const size_t& id) const
{
  CIntEvent eevent(Sim->Dynamics.getEvent(part, Sim->vParticleList[id]));
  if (eevent.getType() != NONE)
    sorter->push(intPart(eevent, eventCount[id]), part.getID());
}

void 
CSNeighbourList::addLocalEvent(const CParticle& part, 
			      const size_t& id) const
{
  if (Sim->Dynamics.getLocals()[id]->isInteraction(part))
    sorter->push(Sim->Dynamics.getLocals()[id]->getEvent(part), part.getID());  
}

void 
CSNeighbourList::addNewEvents(const CParticle& part) const
{  
  //Add the global events
  BOOST_FOREACH(const smrtPlugPtr<CGlobal>& glob, Sim->Dynamics.getGlobals())
    if (glob->isInteraction(part))
      sorter->push(glob->getEvent(part), part.getID());

#ifdef DYNAMO_DEBUG
  if (dynamic_cast<const CGNeighbourList*>
      (Sim->Dynamics.getGlobals()[NBListID].get_ptr())
      == NULL)
    D_throw() << "Not a CGNeighbourList!";
#endif

  //Grab a reference to the neighbour list
  const CGNeighbourList& nblist(*static_cast<const CGNeighbourList*>
				(Sim->Dynamics.getGlobals()[NBListID]
				 .get_ptr()));
  
  //Add the local cell events
  nblist.getParticleLocalNeighbourhood
    (part, boost::bind(&CSNeighbourList::addLocalEvent, this, _1, _2));

  //Add the interaction events
  nblist.getParticleNeighbourhood
    (part, boost::bind(&CSNeighbourList::addInteractionEvent, this, _1, _2));  
}
