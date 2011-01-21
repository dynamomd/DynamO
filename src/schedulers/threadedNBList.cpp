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
  //The operator<<(XML) is virtual but the object is of type
  //CSNeighbourList when it is called
  operator<<(XML);
  I_cout() << "Threaded Variant Loaded with " << _threadPool.getThreadCount()
	   << " threads in the pool";
}

SThreadedNBList::SThreadedNBList(DYNAMO::SimData* const Sim, CSSorter* ns, 
				 size_t threadCount):
  CSNeighbourList(Sim, ns)
{ 
  I_cout() << "Threaded Variant Loaded"; 
  _threadPool.setThreadCount(threadCount);
}

void 
SThreadedNBList::operator<<(const XMLNode& XML)
{
  CSNeighbourList::operator<<(XML);
  _threadPool.setThreadCount(boost::lexical_cast<size_t>
			     (XML.getAttribute("ThreadCount")));
}

void 
SThreadedNBList::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Type") << "ThreadedNeighbourList"
      << xml::attr("ThreadCount") << _threadPool.getThreadCount()
      << xml::tag("Sorter")
      << sorter
      << xml::endtag("Sorter");
}

void 
SThreadedNBList::addEvents(const Particle& part)
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
    (part, magnet::function::MakeDelegate(this, &SThreadedNBList::streamParticles));  

  nblist.getParticleNeighbourhood
    (part, magnet::function::MakeDelegate(this, &SThreadedNBList::addEvents2));  
}

void 
SThreadedNBList::addEventsInit(const Particle& part)
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
    (part, magnet::function::MakeDelegate(this, &CScheduler::addInteractionEventInit));  
}

struct NBlistData {
  void AddNBIDs(const Particle&p1, const size_t& ID) { nbIDs.push_back(ID); }
  
  std::vector<size_t> nbIDs;
};

void 
SThreadedNBList::fullUpdate(const Particle& p1, const Particle& p2)
{
#ifdef DYNAMO_DEBUG
  if (dynamic_cast<const CGNeighbourList*>(Sim->dynamics.getGlobals()[NBListID].get_ptr())
      == NULL)  M_throw() << "Not a CGNeighbourList!";
#endif

  //Now grab a reference to the neighbour list
  const CGNeighbourList& nblist(*static_cast<const CGNeighbourList*>
				(Sim->dynamics.getGlobals()[NBListID].get_ptr()));

  //Now fetch the neighborlist data
  NBlistData nbIDs1, nbIDs2;  
  nblist.getParticleNeighbourhood(p1, magnet::function::MakeDelegate(&nbIDs1, &NBlistData::AddNBIDs));
  nblist.getParticleNeighbourhood(p2, magnet::function::MakeDelegate(&nbIDs2, &NBlistData::AddNBIDs));
  
  //Stream all of the particles up to date
  Sim->dynamics.getLiouvillean().updateParticle(p1);
  Sim->dynamics.getLiouvillean().updateParticle(p2);
  BOOST_FOREACH(const size_t& ID, nbIDs1.nbIDs) Sim->dynamics.getLiouvillean().updateParticle(Sim->particleList[ID]);
  BOOST_FOREACH(const size_t& ID, nbIDs2.nbIDs) Sim->dynamics.getLiouvillean().updateParticle(Sim->particleList[ID]);

  //Both particles events must be invalidated at once
  ++eventCount[p1.getID()];
  ++eventCount[p2.getID()];

  sorter->clearPEL(p1.getID());
  sorter->clearPEL(p2.getID());

  //Add the interaction events, these can churn while try to add the other events
  BOOST_FOREACH(const size_t& ID, nbIDs1.nbIDs) 
    _threadPool.queueTask(magnet::function::Task::makeTask(&SThreadedNBList::threadAddIntEvent, 
						     this, p1, ID, _P1SorterLock));

  BOOST_FOREACH(const size_t& ID, nbIDs2.nbIDs) 
    _threadPool.queueTask(magnet::function::Task::makeTask(&SThreadedNBList::threadAddIntEvent, 
						     this, p2, ID, _P2SorterLock));

  //Add the global events
  BOOST_FOREACH(const magnet::ClonePtr<Global>& glob, Sim->dynamics.getGlobals())
    {
      if (glob->isInteraction(p1))
	_threadPool.queueTask(magnet::function::Task::makeTask(&SThreadedNBList::addGlobal, 
							 this, p1, glob, _P1SorterLock));

      if (glob->isInteraction(p2))
	_threadPool.queueTask(magnet::function::Task::makeTask(&SThreadedNBList::addGlobal, 
							 this, p2, glob, _P2SorterLock));
    }
  
  //Add the local cell events
  nblist.getParticleLocalNeighbourhood
    (p1, magnet::function::MakeDelegate(this, &SThreadedNBList::spawnThreadAddLocalEvent1));

  nblist.getParticleLocalNeighbourhood
    (p2, magnet::function::MakeDelegate(this, &SThreadedNBList::spawnThreadAddLocalEvent2));
    
  //Now wait for the pool
  _threadPool.wait();

  sorter->update(p1.getID());
  sorter->update(p2.getID());
}

void 
SThreadedNBList::addGlobal(const Particle& part, const magnet::ClonePtr<Global>& glob, magnet::thread::Mutex& sorterLock)
{
  GlobalEvent event = glob->getEvent(part);

  magnet::thread::ScopedLock lock1(sorterLock);      
  sorter->push(event, part.getID());
}

void 
SThreadedNBList::fullUpdate(const Particle& part)
{
  invalidateEvents(part);
  addEvents(part);
  sort(part);
}

void 
SThreadedNBList::spawnThreadAddLocalEvent1(const Particle& part, 
					  const size_t& id) 
{
  if (Sim->dynamics.getLocals()[id]->isInteraction(part))
    _threadPool.queueTask(magnet::function::Task::makeTask(&SThreadedNBList::threadAddLocalEvent, 
							   this, part, id, _P1SorterLock));
}

void 
SThreadedNBList::spawnThreadAddLocalEvent2(const Particle& part, 
					  const size_t& id) 
{
  if (Sim->dynamics.getLocals()[id]->isInteraction(part))
    _threadPool.queueTask(magnet::function::Task::makeTask(&SThreadedNBList::threadAddLocalEvent, 
							   this, part, id, _P2SorterLock));
}

void 
SThreadedNBList::threadAddIntEvent(const Particle& part, 
				   const size_t id,
				   magnet::thread::Mutex& sorterLock)
{
  const IntEvent& eevent(Sim->dynamics.getEvent(part, Sim->particleList[id]));
  
  if (eevent.getType() != NONE)
    {
      magnet::thread::ScopedLock lock1(sorterLock);
      sorter->push(intPart(eevent, eventCount[id]), part.getID());
    }
}

void 
SThreadedNBList::threadAddLocalEvent(const Particle& part, 
				     const size_t id,
				     magnet::thread::Mutex& sorterLock)
{
  LocalEvent Event = Sim->dynamics.getLocals()[id]->getEvent(part);

  {
    magnet::thread::ScopedLock lock1(sorterLock);
    sorter->push(Event, part.getID());  
  }
}


void 
SThreadedNBList::threadStreamParticles(const size_t id) const
{
  Sim->dynamics.getLiouvillean().updateParticle(Sim->particleList[id]);
}

void 
SThreadedNBList::streamParticles(const Particle& part, 
				 const size_t& id) const
{
  Sim->dynamics.getLiouvillean().updateParticle(Sim->particleList[id]);
}

void 
SThreadedNBList::addEvents2(const Particle& part, 
			    const size_t& id) const
{
  const IntEvent& eevent(Sim->dynamics.getEvent(part, Sim->particleList[id]));
  
  if (eevent.getType() != NONE)
    sorter->push(intPart(eevent, eventCount[id]), part.getID());
}
