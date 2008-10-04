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

#include "multlist.hpp"
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

void 
CSMultList::stream(const Iflt dt)
{
  eventHeap.stream(dt);
}

const CIntEvent 
CSMultList::earliestIntEvent() const
{
  return Sim->Dynamics.getEvent(Sim->vParticleList[eventHeap.next_ID()], 
				*eventHeap.next_Data().top().p2);
}

const CGlobEvent
CSMultList::earliestGlobEvent() const
{
  CGlobEvent col1, col2;
  BOOST_FOREACH(const smrtPlugPtr<CGlobal>& ptrGlob, Sim->Dynamics.getGlobals())
    if (ptrGlob->isInteraction(Sim->vParticleList[eventHeap.next_ID()]))
      {
	col1 = ptrGlob->getEvent(Sim->vParticleList[eventHeap.next_ID()]);
	if (col1 < col2)
	  col2 = col1;
      }

  return col2;
}

void
CSMultList::initialise()
{
  if (Sim->Dynamics.BCTypeTest<CRLEBC>()
      || Sim->Dynamics.BCTypeTest<CSLEBC>())
    D_throw() << "This scheduler isn't suitable for sheared systems";


  reinitialise(Sim->Dynamics.getLongestInteraction()); 
}

void 
CSMultList::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "MultList";
  CSCells::outputXML(XML);
}

CSMultList::CSMultList(const XMLNode& XML, const DYNAMO::SimData* Sim):
  CSCells(Sim,"MultiList")
{
  
  I_cout() << "Multi List Cellular Algorithmn";
  operator<<(XML);
}

CSMultList::CSMultList(const DYNAMO::SimData* Sim):
  CSCells(Sim,"MultiList")
{ I_cout() << "Multi List Cellular Algorithmn"; }

CSMultList::CSMultList(const DYNAMO::SimData* Sim, const char* Nom):
  CSCells(Sim, Nom)
{}

void 
CSMultList::rescaleTimes(Iflt scale)
{ eventHeap.rescaleTimes(scale); }

void 
CSMultList::update(const CParticle& part)
{
  //Invalidate previous entries
  ++eventCount[part.getID()];
  eventHeap[part.getID()].clear();
  addNewEvents(part);
  eventHeap.update(part.getID());
}

void
CSMultList::reinitialise(Iflt maxdiam)
{
  I_cout() << "Reinitialising on collision " << Sim->lNColl;
  eventHeap.clear();
  eventHeap.resize(Sim->lN);
  eventCount.clear();
  eventCount.resize(Sim->lN, 0);

  //Create the cells
  addCells(maxdiam, false);
  
  //Now initialise the interactions
  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    addNewEvents_init(part);
  
  eventHeap.init();

#ifndef CBT
  I_cout() << "BPQ: Number of lists " << eventHeap.NLists();
  I_cout() << "BPQ: Scale Factor " << eventHeap.scaleFactor();
#endif
}

ENextEvent 
CSMultList::nextEventType() const
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
  
  //First, sort through the interaction events, purging bad events

#ifdef DYNAMO_UpdateCollDebug
  std::cerr << "\nNext eventdt = " << eventHeap.next_dt();
#endif
  while (eventHeap.next_dt() < tmpt)
    {
     switch (eventHeap.next_Data().top().type)
      {
      case INTERACTION:
	if (eventHeap.next_Data().top().collCounter2 != eventCount[eventHeap.next_Data().top().p2->getID()])
	  {
#ifdef DYNAMO_UpdateCollDebug
	    std::cerr << "\nEvent invalid, popping and updating" << eventHeap.next_dt();
#endif
	    eventHeap.next_Data().pop();
	    eventHeap.update(eventHeap.next_ID());
	    break;
	  }

	return Interaction;
      case CELL:
	cellEvent(Sim->vParticleList[eventHeap.next_ID()]);
	break;
      case GLOBAL:
	return Global;
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
CSMultList::addNewEvents_init(const CParticle& part) const
{
  //Add the global event
  if (!Sim->Dynamics.getGlobals().empty())
    eventHeap.push(getGlobEvent(part), part.getID());
  
  //Update the wall collision
  eventHeap.push(Sim->Dynamics.Liouvillean().
		 getSquareCellCollision
		 (part, cells[partCellData[part.getID()].cell].origin, 
		  cellDimension),
		 part.getID());
  
  for (int next = cells[partCellData[part.getID()].cell].list; 
       next != -1; next = partCellData[next].next)
    if (part.getID() != static_cast<unsigned long>(next))
      if (part.getID() < static_cast<unsigned long>(next))
	{
	  CIntEvent eevent(Sim->Dynamics.getEvent(part, Sim->vParticleList[next]));
	  if (eevent.getType() != NONE)
	    eventHeap.push(intPart(eevent, eventCount[next]), part.getID());
	}
  
  //The neighbours particles, same procedure
  BOOST_FOREACH(const int& nb, 
		cells[partCellData[part.getID()].cell].neighbours)
    for (int next = cells[nb].list;
	 next != -1; next = partCellData[next].next)
      if (part.getID() < static_cast<unsigned long>(next))
	{
	  CIntEvent eevent(Sim->Dynamics.getEvent(part, Sim->vParticleList[next]));
	  if (eevent.getType() != NONE)
	    eventHeap.push(intPart(eevent, eventCount[next]), part.getID());
	}
}

void 
CSMultList::addNewEvents(const CParticle& part) const
{  
  //Add the global event
  if (!Sim->Dynamics.getGlobals().empty())
    eventHeap.push(getGlobEvent(part), part.getID());
    
  //Update the wall collision
  eventHeap.push(Sim->Dynamics.Liouvillean().
		 getSquareCellCollision
		 (part, cells[partCellData[part.getID()].cell].origin, 
		  cellDimension),
		 part.getID());

  for (int next = cells[partCellData[part.getID()].cell].list; 
       next != -1; next = partCellData[next].next)
    if (part.getID() != static_cast<unsigned long>(next))
      {
	CIntEvent eevent(Sim->Dynamics.getEvent(part, Sim->vParticleList[next]));
	if (eevent.getType() != NONE)
	  eventHeap.push(intPart(eevent, eventCount[next]), part.getID());
      }

  //The neighbours particles, same procedure
  BOOST_FOREACH(const int& nb, 
		cells[partCellData[part.getID()].cell].neighbours)
    for (int next = cells[nb].list;
	 next != -1; next = partCellData[next].next)
      {
	CIntEvent eevent(Sim->Dynamics.getEvent(part, Sim->vParticleList[next]));
	if (eevent.getType() != NONE)
	  eventHeap.push(intPart(eevent, eventCount[next]), part.getID());
      }
}

void 
CSMultList::cellEvent(const CParticle& part) const 
{  
  //Determine the cell transition direction, its saved
  int cellDirection(eventHeap[part.getID()].top().collCounter2);
  long inPosition;
  int endCell;
  if (std::signbit(part.getVelocity()[cellDirection])) 
    {
      endCell = cells[partCellData[part.getID()].cell].negCells[cellDirection];
      inPosition = cells[cells[endCell].negCells[cellDirection]].coords[cellDirection];
    }
  else
    {
      endCell = cells[partCellData[part.getID()].cell].posCells[cellDirection];
      inPosition = cells[cells[endCell].posCells[cellDirection]].coords[cellDirection];
    }

  //Debug section
#ifdef DYNAMO_WallCollDebug
  {      
    CVector<long> tmp = cells[partCellData[part.getID()].cell].coords;
    CVector<long> tmp2 = cells[endCell].coords;
    
    std::cerr << "\nsysdt " 
	      << (eventHeap.next_dt() + Sim->dSysTime)
	      << "  WALL ID "
	      << part.getID()
	      << "  dt " << (eventHeap.next_dt())
	      << "  from <" 
	      << tmp[0] << "," << tmp[1] << "," << tmp[2]
	      << "> to <" 
	      << tmp2[0] << "," << tmp2[1] << "," << tmp2[2] << ">";
  }
#endif  

  //Remove the cell event from the priority queue (its the top() one)
  eventHeap[part.getID()].pop();  
  //Remove the particle from the cell list
  removeFromCell(part.getID());

  //Particle has just arrived into a new cell, don't do a full update
  //just test it against its new neighbours
  BOOST_FOREACH(const int& nb, cells[endCell].neighbours)
    if (cells[nb].coords[cellDirection] == inPosition)
      for (int next = cells[nb].list; next != -1; next = partCellData[next].next)
	{
	  CIntEvent eevent(Sim->Dynamics.getEvent
			   (part, Sim->vParticleList[next]));
	  if (eevent.getType() != NONE)
	    eventHeap.push(intPart(eevent, eventCount[next]), part.getID());
	}

  //Update Cell Queue
  addToCell(part.getID(), endCell);

  eventHeap.push(Sim->Dynamics.Liouvillean()
		 .getSquareCellCollision(part, cells[endCell].origin, 
					 cellDimension),
		 part.getID());

  eventHeap.update(part.getID());
}

