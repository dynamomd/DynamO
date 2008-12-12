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
#include "../simulation/particle.hpp"
#include "../dynamics/dynamics.hpp"
#include "../dynamics/liouvillean/liouvillean.hpp"
#include "../base/is_simdata.hpp"
#include "../base/is_base.hpp"
#include "../dynamics/systems/system.hpp"
#include <cmath> //for huge val
#include "../extcode/xmlParser.h"
#include "../dynamics/globals/global.hpp"
#include "../dynamics/globals/globEvent.hpp"
#include "../dynamics/globals/neighbourList.hpp"
#include "../dynamics/locals/local.hpp"
#include "../dynamics/locals/localEvent.hpp"
#include <boost/bind.hpp>

void 
CSNeighbourList::operator<<(const XMLNode& XML)
{
  sorter.set_ptr(CSSorter::getClass(XML.getChildNode("Sorter"), Sim));
}

void
CSNeighbourList::initialise()
{
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
  //The plus one is because system events are stored in the last heap;
  sorter->resize(Sim->lN+1);
  eventCount.clear();
  eventCount.resize(Sim->lN+1, 0);

  //Now initialise the interactions
  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    addEvents(part);
  
  sorter->init();

  rebuildSystemEvents();
  
  //Register the new neighbour function with the cellular tracker
  if (!cellChange)
    cellChange =
      static_cast<const CGNeighbourList&>
      (*(Sim->Dynamics.getGlobals()[NBListID]))
      .ConnectSigNewNeighbourNotify
      (&CSNeighbourList::addInteractionEvent, this);

  if (!cellChangeLocal)
    cellChangeLocal =
      static_cast<const CGNeighbourList&>
      (*(Sim->Dynamics.getGlobals()[NBListID]))
      .ConnectSigNewLocalNotify
      (&CSNeighbourList::addLocalEvent, this);
  
  if (!reinit)
    reinit = 
      static_cast<const CGNeighbourList&>
      (*(Sim->Dynamics.getGlobals()[NBListID]))
      .ConnectSigReInitNotify
      (&CSNeighbourList::initialise, this);
}

void 
CSNeighbourList::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "NeighbourList"
      << xmlw::tag("Sorter")
      << sorter
      << xmlw::endtag("Sorter");
}

CSNeighbourList::CSNeighbourList(const XMLNode& XML, 
				 const DYNAMO::SimData* Sim):
  CScheduler(Sim,"NeighbourListScheduler", NULL),
  cellChange(0),
  cellChangeLocal(0),
  reinit(0)
{ 
  I_cout() << "Neighbour List Scheduler Algorithmn Loaded";
  operator<<(XML);
}

CSNeighbourList::CSNeighbourList(const CSNeighbourList& nb):
  CScheduler(nb),
  NBListID(nb.NBListID),
  cellChange(0),
  cellChangeLocal(0),
  reinit(0)
{}

CSNeighbourList::CSNeighbourList(const DYNAMO::SimData* Sim, CSSorter* ns):
  CScheduler(Sim,"NeighbourListScheduler", ns),
  cellChange(0),
  cellChangeLocal(0),
  reinit(0)
{ I_cout() << "Neighbour List Scheduler Algorithmn Loaded"; }

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
CSNeighbourList::addEvents(const CParticle& part)
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
    (part, CGNeighbourList::getNBDelegate(&CSNeighbourList::addLocalEvent, this));

  //Add the interaction events
  nblist.getParticleNeighbourhood
    (part, CGNeighbourList::getNBDelegate(&CSNeighbourList::addInteractionEvent, this));  
}
