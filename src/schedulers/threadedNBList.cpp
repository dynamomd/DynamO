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

#include "threadedNBList.hpp"
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
#include <boost/progress.hpp>

SThreadedNBList::SThreadedNBList(const XMLNode& XML, 
				 DYNAMO::SimData* const Sim):
  CSNeighbourList(XML, Sim)
{ 
  I_cout() << "Threaded Variant Loaded";
  _threadPool.setThreadCount(4);
}

SThreadedNBList::SThreadedNBList(DYNAMO::SimData* const Sim, CSSorter* ns):
  CSNeighbourList(Sim, ns)
{ 
  I_cout() << "Threaded Variant Loaded"; 
  _threadPool.setThreadCount(4);
}


void 
SThreadedNBList::operator<<(const XMLNode& XML)
{
  CSNeighbourList::operator<<(XML);
}

void 
SThreadedNBList::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "ThreadedNeighbourList"
      << xmlw::tag("Sorter")
      << sorter
      << xmlw::endtag("Sorter");
}

void 
SThreadedNBList::addEvents(const CParticle& part)
{
  Sim->dynamics.getLiouvillean().updateParticle(part);
  
  //Add the global events
  BOOST_FOREACH(const smrtPlugPtr<CGlobal>& glob, Sim->dynamics.getGlobals())
    if (glob->isInteraction(part))
      sorter->push(glob->getEvent(part), part.getID());
  
#ifdef DYNAMO_DEBUG
  if (dynamic_cast<const CGNeighbourList*>
      (Sim->dynamics.getGlobals()[NBListID].get_ptr())
      == NULL)
    D_throw() << "Not a CGNeighbourList!";
#endif

  //Grab a reference to the neighbour list
  const CGNeighbourList& nblist(*static_cast<const CGNeighbourList*>
				(Sim->dynamics.getGlobals()[NBListID]
				 .get_ptr()));
  
  //Add the local cell events
  nblist.getParticleLocalNeighbourhood
    (part, fastdelegate::MakeDelegate(this, &CScheduler::addLocalEvent));

  //Add the interaction events
  nblist.getParticleNeighbourhood
    (part, fastdelegate::MakeDelegate(this, &SThreadedNBList::streamParticles));  

  nblist.getParticleNeighbourhood
    (part, fastdelegate::MakeDelegate(this, &SThreadedNBList::addEvents2));  
}

void 
SThreadedNBList::streamParticles(const CParticle& part, 
				 const size_t& id) const
{
  Sim->dynamics.getLiouvillean().updateParticle(Sim->vParticleList[id]);
}

void 
SThreadedNBList::addEvents2(const CParticle& part, 
			    const size_t& id) const
{
  const CIntEvent& eevent(Sim->dynamics.getEvent(part, Sim->vParticleList[id]));
  
  if (eevent.getType() != NONE)
    sorter->push(intPart(eevent, eventCount[id]), part.getID());
}


void 
SThreadedNBList::addEventsInit(const CParticle& part)
{  
  Sim->dynamics.getLiouvillean().updateParticle(part);

  //Add the global events
  BOOST_FOREACH(const smrtPlugPtr<CGlobal>& glob, Sim->dynamics.getGlobals())
    if (glob->isInteraction(part))
      sorter->push(glob->getEvent(part), part.getID());
  
#ifdef DYNAMO_DEBUG
  if (dynamic_cast<const CGNeighbourList*>
      (Sim->dynamics.getGlobals()[NBListID].get_ptr())
      == NULL)
    D_throw() << "Not a CGNeighbourList!";
#endif

  //Grab a reference to the neighbour list
  const CGNeighbourList& nblist(*static_cast<const CGNeighbourList*>
				(Sim->dynamics.getGlobals()[NBListID]
				 .get_ptr()));
  
  //Add the local cell events
  nblist.getParticleLocalNeighbourhood
    (part, fastdelegate::MakeDelegate(this, &CScheduler::addLocalEvent));

  //Add the interaction events
  nblist.getParticleNeighbourhood
    (part, fastdelegate::MakeDelegate(this, &CScheduler::addInteractionEventInit));  
}

void 
SThreadedNBList::fullUpdate(const CParticle& p1, const CParticle& p2)
{
  //Both must be invalidated at once to reduce the number of invalid
  //events in the queue
  ++eventCount[p1.getID()];
  ++eventCount[p2.getID()];

  sorter->clearPEL(p1.getID());
  Sim->dynamics.getLiouvillean().updateParticle(p1);
  sorter->clearPEL(p2.getID());
  Sim->dynamics.getLiouvillean().updateParticle(p2);

  //Add the global events
  BOOST_FOREACH(const smrtPlugPtr<CGlobal>& glob, Sim->dynamics.getGlobals())
    if (glob->isInteraction(p1))
      _threadPool.invoke(boost::bind(&SThreadedNBList::addGlobal, this, boost::ref(p1), boost::ref(glob)));

  BOOST_FOREACH(const smrtPlugPtr<CGlobal>& glob, Sim->dynamics.getGlobals())
    if (glob->isInteraction(p2))
      _threadPool.invoke(boost::bind(&SThreadedNBList::addGlobal, this, boost::ref(p2), boost::ref(glob)));
  
  _threadPool.wait();

#ifdef DYNAMO_DEBUG
  if (dynamic_cast<const CGNeighbourList*>
      (Sim->dynamics.getGlobals()[NBListID].get_ptr())
      == NULL)
    D_throw() << "Not a CGNeighbourList!";
#endif

  //Grab a reference to the neighbour list
  const CGNeighbourList& nblist(*static_cast<const CGNeighbourList*>
				(Sim->dynamics.getGlobals()[NBListID]
				 .get_ptr()));
  
  //Add the local cell events
  nblist.getParticleLocalNeighbourhood
    (p1, fastdelegate::MakeDelegate(this, &CScheduler::addLocalEvent));
  nblist.getParticleLocalNeighbourhood
    (p2, fastdelegate::MakeDelegate(this, &CScheduler::addLocalEvent));

  //Add the interaction events
  nblist.getParticleNeighbourhood
    (p1, fastdelegate::MakeDelegate(this, &SThreadedNBList::streamParticles));  

  nblist.getParticleNeighbourhood
    (p2, fastdelegate::MakeDelegate(this, &SThreadedNBList::streamParticles));  

  nblist.getParticleNeighbourhood
    (p1, fastdelegate::MakeDelegate(this, &SThreadedNBList::addEvents2));  

  nblist.getParticleNeighbourhood
    (p2, fastdelegate::MakeDelegate(this, &SThreadedNBList::addEvents2));  

  sorter->update(p1.getID());
  sorter->update(p2.getID());
}

void 
SThreadedNBList::addGlobal(const CParticle& part, const smrtPlugPtr<CGlobal>& glob)
{
  CGlobEvent event = glob->getEvent(part);

  boost::mutex::scoped_lock lock1(_sorterLock);      
  sorter->push(event, part.getID());
}

void 
SThreadedNBList::fullUpdate(const CParticle& part)
{
  invalidateEvents(part);
  addEvents(part);
  sort(part);
}
