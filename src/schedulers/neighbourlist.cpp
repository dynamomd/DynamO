/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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
#include "../dynamics/globals/global.hpp"
#include "../dynamics/globals/globEvent.hpp"
#include "../dynamics/globals/neighbourList.hpp"
#include "../dynamics/locals/local.hpp"
#include "../dynamics/locals/localEvent.hpp"
#include <magnet/xmlreader.hpp>
#include <boost/bind.hpp>
#include <boost/progress.hpp>
#include <cmath> //for huge val

void 
CSNeighbourList::operator<<(const magnet::xml::Node& XML)
{
  sorter.set_ptr(CSSorter::getClass(XML.getNode("Sorter"), Sim));
}

void
CSNeighbourList::initialise()
{
  try {
    NBListID = Sim->dynamics.getGlobal("SchedulerNBList")->getID();
  }
  catch(std::exception& cxp)
    {
      M_throw() << "Failed while finding the neighbour list global.\n"
		<< "You must have a neighbour list enabled for this\n"
		<< "scheduler called SchedulerNBList.\n"
		<< cxp.what();
    }
  
  if (dynamic_cast<const CGNeighbourList*>
      (Sim->dynamics.getGlobals()[NBListID].get_ptr())
      == NULL)
    M_throw() << "The Global named SchedulerNBList is not a neighbour list!";

  static_cast<CGNeighbourList&>
    (*Sim->dynamics.getGlobals()[NBListID].get_ptr())
    .markAsUsedInScheduler();

  I_cout() << "Building all events on collision " << Sim->eventCount;
  std::cout.flush();

  sorter->clear();
  //The plus one is because system events are stored in the last heap;
  sorter->resize(Sim->N+1);
  eventCount.clear();
  eventCount.resize(Sim->N+1, 0);

  //Now initialise the interactions
  {
    boost::progress_display prog(Sim->N);
 
    BOOST_FOREACH(const Particle& part, Sim->particleList)
      {
	addEventsInit(part);
	++prog;
      }
  }
  
  sorter->init();

  rebuildSystemEvents();
}

void 
CSNeighbourList::rebuildList()
{ 
#ifdef DYNAMO_DEBUG
  initialise();
#else
  sorter->clear();
  //The plus one is because system events are stored in the last heap;
  sorter->resize(Sim->N+1);
  eventCount.clear();
  eventCount.resize(Sim->N+1, 0);

  BOOST_FOREACH(const Particle& part, Sim->particleList)
    addEventsInit(part);
  
  sorter->rebuild();
  
  rebuildSystemEvents();
#endif
}

void 
CSNeighbourList::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Type") << "NeighbourList"
      << xml::tag("Sorter")
      << sorter
      << xml::endtag("Sorter");
}

CSNeighbourList::CSNeighbourList(const magnet::xml::Node& XML, 
				 DYNAMO::SimData* const Sim):
  CScheduler(Sim,"NeighbourListScheduler", NULL)
{ 
  I_cout() << "Neighbour List Scheduler Algorithmn Loaded";
  operator<<(XML);
}

CSNeighbourList::CSNeighbourList(DYNAMO::SimData* const Sim, CSSorter* ns):
  CScheduler(Sim,"NeighbourListScheduler", ns)
{ I_cout() << "Neighbour List Scheduler Algorithmn Loaded"; }

void 
CSNeighbourList::addEvents(const Particle& part)
{
  Sim->dynamics.getLiouvillean().updateParticle(part);
  
  //Add the global events
  BOOST_FOREACH(const magnet::ClonePtr<Global>& glob, Sim->dynamics.getGlobals())
    if (glob->isInteraction(part))
      sorter->push(glob->getEvent(part), part.getID());
  
#ifdef DYNAMO_DEBUG
  if (dynamic_cast<const CGNeighbourList*>
      (Sim->dynamics.getGlobals()[NBListID].get_ptr())
      == NULL)
    M_throw() << "Not a CGNeighbourList!";
#endif

  //Grab a reference to the neighbour list
  const CGNeighbourList& nblist(*static_cast<const CGNeighbourList*>
				(Sim->dynamics.getGlobals()[NBListID]
				 .get_ptr()));
  
  //Add the local cell events
  nblist.getParticleLocalNeighbourhood
    (part, magnet::function::MakeDelegate(this, &CScheduler::addLocalEvent));

  //Add the interaction events
  nblist.getParticleNeighbourhood
    (part, magnet::function::MakeDelegate(this, &CScheduler::addInteractionEvent));  
}

void 
CSNeighbourList::addEventsInit(const Particle& part)
{  
  Sim->dynamics.getLiouvillean().updateParticle(part);

  //Add the global events
  BOOST_FOREACH(const magnet::ClonePtr<Global>& glob, Sim->dynamics.getGlobals())
    if (glob->isInteraction(part))
      sorter->push(glob->getEvent(part), part.getID());
  
#ifdef DYNAMO_DEBUG
  if (dynamic_cast<const CGNeighbourList*>
      (Sim->dynamics.getGlobals()[NBListID].get_ptr())
      == NULL)
    M_throw() << "Not a CGNeighbourList!";
#endif

  //Grab a reference to the neighbour list
  const CGNeighbourList& nblist(*static_cast<const CGNeighbourList*>
				(Sim->dynamics.getGlobals()[NBListID]
				 .get_ptr()));
  
  //Add the local cell events
  nblist.getParticleLocalNeighbourhood
    (part, magnet::function::MakeDelegate
     (this, &CScheduler::addLocalEvent));

  //Add the interaction events
  nblist.getParticleNeighbourhood
    (part, magnet::function::MakeDelegate
     (this, &CScheduler::addInteractionEventInit));  
}
