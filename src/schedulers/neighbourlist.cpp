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
#include "../dynamics/BC/BC.hpp"
#include "../dynamics/BC/LEBC.hpp"
#include "../base/is_simdata.hpp"
#include "../base/is_base.hpp"
#include "../dynamics/globals/globEvent.hpp"
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
      = static_cast<const CGNeighbourList&>
      (*(Sim->Dynamics.getGlobals()[NBListID]))
      .registerCellTransitionNewNeighbourCallBack
      (boost::bind(&CSNeighbourList::addInteractionEvent, this, _1, _2));

  if (!cellChangeLocal.connected())
    cellChangeLocal 
      = static_cast<const CGNeighbourList&>
      (*(Sim->Dynamics.getGlobals()[NBListID]))
      .registerCellTransitionNewLocalCallBack
      (boost::bind(&CSNeighbourList::addLocalEvent, this, _1, _2));
  
  if (!reinit.connected())
    reinit 
      = static_cast<const CGNeighbourList&>
      (*(Sim->Dynamics.getGlobals()[NBListID]))
      .registerReInitNotify
      (boost::bind(&CSNeighbourList::initialise, this));
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
  CScheduler(Sim,"NeighbourListScheduler", NULL)
{ 
  I_cout() << "Neighbour List Scheduler Algorithmn Loaded";
  operator<<(XML);
}

CSNeighbourList::CSNeighbourList(const DYNAMO::SimData* Sim, CSSorter* ns):
  CScheduler(Sim,"NeighbourListScheduler", ns)
{ I_cout() << "Neighbour List Scheduler Algorithmn Loaded"; }

void 
CSNeighbourList::update(const CParticle& part)
{
  //Invalidate previous entries
  ++eventCount[part.getID()];
  (*sorter)[part.getID()].clear();
  addNewEvents(part);
  sorter->update(part.getID());
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
