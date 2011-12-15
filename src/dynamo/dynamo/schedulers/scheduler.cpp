/*  dynamo:- Event driven molecular dynamics simulator 
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

#include <dynamo/schedulers/include.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/dynamics/interactions/intEvent.hpp>
#include <dynamo/dynamics/globals/global.hpp>
#include <dynamo/dynamics/globals/globEvent.hpp>
#include <dynamo/dynamics/locals/local.hpp>
#include <dynamo/dynamics/systems/system.hpp>
#include <dynamo/dynamics/liouvillean/liouvillean.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <dynamo/dynamics/units/units.hpp>

#ifdef DYNAMO_DEBUG
#include <dynamo/dynamics/globals/neighbourList.hpp>
#include <dynamo/dynamics/NparticleEventData.hpp>
#endif

#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

namespace dynamo {
  Scheduler::Scheduler(dynamo::SimData* const tmp, const char * aName,
			 FEL* nS):
    SimBase(tmp, aName),
    sorter(nS),
    _interactionRejectionCounter(0),
    _localRejectionCounter(0)
  {}

  Scheduler::~Scheduler() {}

  void 
  Scheduler::operator<<(const magnet::xml::Node& XML)
  {
    sorter = FEL::getClass(XML.getNode("Sorter"), Sim);
  }

  void
  Scheduler::initialise()
  {
    dout << "Building all events on collision " << Sim->eventCount << std::endl;
    rebuildList();
  }

  void
  Scheduler::rebuildList()
  {
    sorter->clear();
    //The plus one is because system events are stored in the last heap;
    sorter->resize(Sim->N+1);
    eventCount.clear();
    eventCount.resize(Sim->N+1, 0);

    BOOST_FOREACH(const Particle& part, Sim->particleList)
      addEvents(part);
  
    sorter->init();

    rebuildSystemEvents();
  }


  void 
  Scheduler::addEvents(const Particle& part)
  {  
    Sim->dynamics.getLiouvillean().updateParticle(part);

    //Add the global events
    BOOST_FOREACH(const shared_ptr<Global>& glob, Sim->dynamics.getGlobals())
      if (glob->isInteraction(part))
	sorter->push(glob->getEvent(part), part.getID());
  
    //Add the local cell events
    getLocalNeighbourhood
      (part, magnet::function::MakeDelegate(this, &Scheduler::addLocalEvent));

    //Add the interaction events
    getParticleNeighbourhood
      (part, magnet::function::MakeDelegate(this, &Scheduler::addInteractionEvent));
  }

  shared_ptr<Scheduler>
  Scheduler::getClass(const magnet::xml::Node& XML, dynamo::SimData* const Sim)
  {
    if (!strcmp(XML.getAttribute("Type"),"NeighbourList"))
      return shared_ptr<Scheduler>(new SNeighbourList(XML, Sim));
    else if (!strcmp(XML.getAttribute("Type"),"Dumb"))
      return shared_ptr<Scheduler>(new SDumb(XML, Sim));
    else if (!strcmp(XML.getAttribute("Type"),"SystemOnly"))
      return shared_ptr<Scheduler>(new SSystemOnly(XML, Sim));
    else if (!strcmp(XML.getAttribute("Type"),"Complex"))
      return shared_ptr<Scheduler>(new SComplex(XML, Sim));
    else 
      M_throw() << XML.getAttribute("Type")
		<< ", Unknown type of Scheduler encountered";
  }

  magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream& XML, 
				     const Scheduler& g)
  {
    g.outputXML(XML);
    return XML;
  }

  void 
  Scheduler::rebuildSystemEvents() const
  {
    sorter->clearPEL(Sim->N);

    BOOST_FOREACH(const shared_ptr<System>& sysptr, 
		  Sim->dynamics.getSystemEvents())
      sorter->push(Event(sysptr->getdt(), SYSTEM, sysptr->getID(), 0), Sim->N);

    sorter->update(Sim->N);
  }

  void 
  Scheduler::popNextEvent()
  {
    sorter->popNextPELEvent(sorter->next_ID());
  }

  void 
  Scheduler::pushEvent(const Particle& part,
		       const Event& newevent)
  {
    sorter->push(newevent, part.getID());
  }

  void 
  Scheduler::sort(const Particle& part)
  {
    sorter->update(part.getID());
  }

  void 
  Scheduler::invalidateEvents(const Particle& part)
  {
    //Invalidate previous entries
    ++eventCount[part.getID()];
    sorter->clearPEL(part.getID());
  }

  void
  Scheduler::runNextEvent()
  {
    sorter->sort();

#ifdef DYNAMO_DEBUG
    if (sorter->nextPELEmpty())
      M_throw() << "Next particle list is empty but top of list!";
#endif

    lazyDeletionCleanup();

    if (boost::math::isnan(sorter->next_dt()))
      M_throw() << "Next event time is NaN"
		<< "\nTime to event "
		<< sorter->next_dt()
		<< "\nEvent Type = " 
		<< sorter->next_type()
		<< "\nOwner Particle = " << sorter->next_ID()
		<< "\nID2 = " << sorter->next_p2();
    
    if (sorter->next_dt() == HUGE_VAL)
      {
	derr << "Next event time is Inf! (Queue has run out of events!)\n"
	     << "Shutting simulation down..."
	     << "\nEvent details, Type = " 
	     << sorter->next_type()
	     << "\nOwner Particle = " << sorter->next_ID()
	     << "\nID2 = " << sorter->next_p2()
	     << std::endl;
	Sim->endEventCount = Sim->eventCount;
	return;
      }

    ////////////////////////////////////////////////////////////////////
    // We can't perform such strict testing as commented out
    // below. Sometimes negative event times occur, usually at the start
    // of a simulation when particles are initialized just on the edge
    // of a cell, or if we have a system event which is "triggered" and
    // sets its own event time to 0. These must be tolerated and we must
    // trust in the determinism of the dynamics and the precision of the
    // calculations to minimise any effects. Generally, systems
    // shouldn't crash because of negative event times that were not
    // caused by a physically incorrect initial configuration
    ////////////////////////////////////////////////////////////////////
    //  if (sorter->next_dt() < 0)
    //    M_throw() << "Next event time is less than 0"
    //	      << "\nTime to event "
    //	      << sorter->next_dt()
    //	      << "\nEvent Type = " 
    //	      << sorter->next_type()
    //      	      << "\nOwner Particle = " << sorter->next_ID()
    //	      << "\nID2 = " << sorter->next_p2();
  
    /*! This is our dimensionless parameter which we need to correct a
      edge case for the collision testing. If an event is scheduled to
      occur its collision time is always double checked before it is
      executed. If two events are close together in time, the earliest
      might be popped off the queue, retested and then appear to occur
      later than the next event. In this case the original event is
      discarded and the new version is reinserted into the event
      queue. However, a rounding error might then cause the new event to
      appear earlier than the second event and we're back where we
      started. Basically, if "rejectionLimit" rejections occur in a row we
      just accept the next event in the queue. This breaks these loops and
      allows the simulation to continue. 

      With this method the system is guarranteed to maintain the correct
      event sequence to within machine precision. The queue can even
      handle negative time events provided the dynamics allow it.
    */
    const size_t rejectionLimit = 10;

    switch (sorter->next_type())
      {
      case INTERACTION:
	{
	  const Particle& p1(Sim->particleList[sorter->next_ID()]);
	  const Particle& p2(Sim->particleList[sorter->next_p2()]);

	  //Ready the next event in the FEL
	  sorter->popNextEvent();
	  sorter->update(sorter->next_ID());
	  sorter->sort();	
	  lazyDeletionCleanup();

	  //Now recalculate the FEL event
	  Sim->dynamics.getLiouvillean().updateParticlePair(p1, p2);       
	  IntEvent Event(Sim->dynamics.getEvent(p1, p2));
	
#ifdef DYNAMO_DEBUG
	  if (sorter->nextPELEmpty())
	    M_throw() << "The next PEL is empty, cannot perform the comparison to see if this event is out of sequence";
#endif

	  if ((Event.getType() == NONE)
	      || ((Event.getdt() > sorter->next_dt()) 
		  && (++_interactionRejectionCounter < rejectionLimit)))
	    {
#ifdef DYNAMO_DEBUG
	      derr << "Event " << Sim->eventCount << ":" << Event.getType()
		   << ",dt=" << Event.getdt() << ">nextdt=" << sorter->next_dt()
		   << ",p1=" << p1.getID() << ",p2=" << p2.getID() << std::endl;
#endif		
	      this->fullUpdate(p1, p2);
	      return;
	    }
	
#ifdef DYNAMO_DEBUG
	  if (Event.getdt() < 0)
	    derr << "Negative time " << Event.getdt() << std::endl;
#endif

	  //Reset the rejection watchdog, we will run an interaction event now
	  _interactionRejectionCounter = 0;
		
#ifdef DYNAMO_DEBUG

	  if (boost::math::isnan(Event.getdt()))
	    M_throw() << "A NAN Interaction collision time has been found"
		      << Event.stringData(Sim);
	
	  if (Event.getdt() == HUGE_VAL)
	    M_throw() << "An infinite Interaction (not marked as NONE) collision time has been found\n"
		      << Event.stringData(Sim);
#endif
	
	  //Debug section
#ifdef DYNAMO_CollDebug
	  std::cerr << "\nsysdt " << (Event.getdt() + Sim->dSysTime) / Sim->dynamics.units().unitTime()
		    << "  ID1 " << ((p2.getID() < p2.getID()) ? p1.getID() : p2.getID())
		    << "  ID2 " << ((p2.getID() < p2.getID()) ? p2.getID() : p1.getID())
		    << "  dt " << Event.getdt()
		    << "  Type " << Event.getType();
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

	  //We also don't recheck Global events! (Check, some events might rely on this behavior)
	  Sim->dynamics.getGlobals()[sorter->next_p2()]
	    ->runEvent(Sim->particleList[sorter->next_ID()], sorter->next_dt());       	
	  break;	           
	}
      case LOCAL:
	{
	  const Particle& part(Sim->particleList[sorter->next_ID()]);

	  //Copy the FEL event
	  size_t localID = sorter->next_p2();

	  //Ready the next event in the FEL
	  sorter->popNextEvent();
	  sorter->update(sorter->next_ID());
	  sorter->sort();
	  lazyDeletionCleanup();

	  Sim->dynamics.getLiouvillean().updateParticle(part);
	  LocalEvent iEvent(Sim->dynamics.getLocals()[localID]->getEvent(part));

	  if (iEvent.getType() == NONE)
	    {
#ifdef DYNAMO_DEBUG
	      derr << "Local event found not to occur [" << part.getID()
		   << "] (possible glancing/tenuous event canceled due to numerical error)" << std::endl;
#endif		
	      this->fullUpdate(part);
	      return;
	    }
	
	  double next_dt = sorter->next_dt();

	  if ((iEvent.getdt() > next_dt) && (++_localRejectionCounter < rejectionLimit))
	    {
#ifdef DYNAMO_DEBUG 
	      derr << "Recalculated LOCAL event time is greater than the next event time, recalculating" << std::endl;
#endif
	      this->fullUpdate(part);
	      return;
	    }

	  _localRejectionCounter = 0;

#ifdef DYNAMO_DEBUG 
	  if (boost::math::isnan(iEvent.getdt()))
	    M_throw() << "A NAN Global collision time has been found\n"
		      << iEvent.stringData(Sim);
	
	  if (iEvent.getdt() == HUGE_VAL)
	    M_throw() << "An infinite (not marked as NONE) Global collision time has been found\n"
		      << iEvent.stringData(Sim);
#endif
	
	  Sim->dSysTime += iEvent.getdt();
	
	  stream(iEvent.getdt());
	
	  //dynamics must be updated first
	  Sim->dynamics.stream(iEvent.getdt());
	
	  iEvent.addTime(Sim->freestreamAcc);
	  Sim->freestreamAcc = 0;

	  Sim->dynamics.getLocals()[localID]->runEvent(part, iEvent);	  
	  break;
	}
      case SYSTEM:
	{
	  Sim->dynamics.getSystemEvents()[sorter->next_p2()]
	    ->runEvent();
	  //This saves the system events rebuilding themselves
	  rebuildSystemEvents();
	  break;
	}
      case VIRTUAL:
	{
	  //This is a special type which requires that the system is
	  //moved forward to the current time and the events for this
	  //particle recalculated.

	  double dt = sorter->next_dt();
	  size_t ID = sorter->next_ID();
	  Sim->dSysTime += dt;
	  Sim->ptrScheduler->stream(dt);
	  Sim->dynamics.stream(dt);
	  Sim->freestreamAcc += dt;
	  this->fullUpdate(Sim->particleList[ID]);
	  break;
	}
      case NONE:
	{
	  M_throw() << "A NONE event has reached the top of the queue."
	    "\nThe simulation has run out of events! Aborting!";
	}
      default:
	M_throw() << "Unhandled event type requested to be run\n"
		  << "Type is " << sorter->next_type();
      }
  }

  void 
  Scheduler::addInteractionEvent(const Particle& part, 
				 const size_t& id) const
  {
    if (part.getID() == id) return;
    const Particle& part2(Sim->particleList[id]);

    Sim->dynamics.getLiouvillean().updateParticle(part2);

    const IntEvent& eevent(Sim->dynamics.getEvent(part, part2));

    if (eevent.getType() != NONE)
      sorter->push(Event(eevent, eventCount[id]), part.getID());
  }

  void 
  Scheduler::addInteractionEventInit(const Particle& part, 
				      const size_t& id) const
  {
    if (part.getID() == id) return;

    //We'll be smart about memory and try to add events evenly on
    //initialisation to all particles
    //
    //This is achieved by only allowing one particle to store the
    //event. But we can't just use sorting (e.g., part.getID() < id) to
    //discriminate which particle gets the event as, on regular
    //lattices, particle 0 will get all of the events for all the
    //particles in its neighborhood!
    //
    //We can mix this up a little, by also testing for odd and evenness
    //and using this to switch which particles are chosen

    size_t val = (part.getID() % 2) + 2 * (id % 2);
    switch (val)//part-id
      {
      case 0: //even-even (accept half)
	if (part.getID() > id) return;      
	break;
      case 1: //odd-even (accept)
	break;
      case 2: //even-odd (reject)
	return;
      case 3: //odd-odd (accept half)
	if (part.getID() < id) return;      
      }

    addInteractionEvent(part, id);
  }

  void 
  Scheduler::addLocalEvent(const Particle& part, 
			   const size_t& id) const
  {
    if (Sim->dynamics.getLocals()[id]->isInteraction(part))
      sorter->push(Sim->dynamics.getLocals()[id]->getEvent(part), part.getID());  
  }

  void 
  Scheduler::lazyDeletionCleanup()
  {
    while ((sorter->next_type() == INTERACTION)
	   && (sorter->next_collCounter2()
	       != eventCount[sorter->next_p2()]))
      {
	//Not valid, update the list
	sorter->popNextEvent();
	sorter->update(sorter->next_ID());
	sorter->sort();
      
#ifdef DYNAMO_DEBUG
	if (sorter->nextPELEmpty())
	  M_throw() << "Next particle list is empty but top of list!";
#endif
      }
  }
}
