/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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
#include <dynamo/interactions/intEvent.hpp>
#include <dynamo/globals/global.hpp>
#include <dynamo/globals/globEvent.hpp>
#include <dynamo/locals/local.hpp>
#include <dynamo/systems/system.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/units/units.hpp>

#ifdef DYNAMO_DEBUG
#include <dynamo/globals/neighbourList.hpp>
#include <dynamo/NparticleEventData.hpp>
#endif

#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

namespace dynamo {
  Scheduler::Scheduler(dynamo::Simulation* const tmp, const char * aName,
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
    //Now, the scheduler is used to test the state of the system.
    dout << "Checking the simulation configuration for any errors" << std::endl;
    size_t warnings(0);

    BOOST_FOREACH(const shared_ptr<Interaction>& interaction_ptr, Sim->interactions)
      warnings += interaction_ptr->validateState(warnings < 101, 101 - warnings);
    
    for (size_t id1(0); id1 < Sim->particles.size(); ++id1)
      {
	std::auto_ptr<IDRange> ids(getParticleNeighbours(Sim->particles[id1]));
	BOOST_FOREACH(const size_t id2, *ids)
	  if (id2 > id1)
	    if (Sim->getInteraction(Sim->particles[id1], Sim->particles[id2])
		->validateState(Sim->particles[id1], Sim->particles[id2], (warnings < 101)))
	      ++warnings;
      }
    
    BOOST_FOREACH(const Particle& part, Sim->particles)
      BOOST_FOREACH(const shared_ptr<Local>& lcl, Sim->locals)
      if (lcl->isInteraction(part))
	if (lcl->validateState(part, (warnings < 101)))
	  ++warnings;
    
    if (warnings > 100)
      derr << "Over 100 warnings of invalid states, further output was suppressed (total of " << warnings << " warnings detected)" << std::endl;

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

    BOOST_FOREACH(Particle& part, Sim->particles)
      addEvents(part);
  
    sorter->init();

    rebuildSystemEvents();
  }


  void 
  Scheduler::addEvents(Particle& part)
  {  
    Sim->dynamics->updateParticle(part);

    //Add the global events
    BOOST_FOREACH(const shared_ptr<Global>& glob, Sim->globals)
      if (glob->isInteraction(part))
	sorter->push(glob->getEvent(part), part.getID());
  
    //Add the local cell events
    std::auto_ptr<IDRange> ids(getParticleLocals(part));
    
    BOOST_FOREACH(const size_t id2, *ids)
      addLocalEvent(part, id2);

    //Now add the interaction events
    ids = getParticleNeighbours(part);
    BOOST_FOREACH(const size_t id2, *ids)
      addInteractionEvent(part, id2);
  }

  shared_ptr<Scheduler>
  Scheduler::getClass(const magnet::xml::Node& XML, dynamo::Simulation* const Sim)
  {
    if (!XML.getAttribute("Type").getValue().compare("NeighbourList"))
      return shared_ptr<Scheduler>(new SNeighbourList(XML, Sim));
    else if (!XML.getAttribute("Type").getValue().compare("Dumb"))
      return shared_ptr<Scheduler>(new SDumb(XML, Sim));
    else if (!XML.getAttribute("Type").getValue().compare("SystemOnly"))
      return shared_ptr<Scheduler>(new SSystemOnly(XML, Sim));
    else if (!XML.getAttribute("Type").getValue().compare("Complex"))
      return shared_ptr<Scheduler>(new SComplex(XML, Sim));
    else 
      M_throw() << XML.getAttribute("Type").getValue()
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
		  Sim->systems)
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
	  Particle& p1(Sim->particles[sorter->next_ID()]);
	  Particle& p2(Sim->particles[sorter->next_p2()]);

	  //Ready the next event in the FEL
	  sorter->popNextEvent();
	  sorter->update(sorter->next_ID());
	  sorter->sort();	
	  lazyDeletionCleanup();

	  //Now recalculate the FEL event
	  Sim->dynamics->updateParticlePair(p1, p2);       
	  IntEvent Event(Sim->getEvent(p1, p2));
	
#ifdef DYNAMO_DEBUG
	  if (sorter->nextPELEmpty())
	    M_throw() << "The next PEL is empty, cannot perform the comparison to see if this event is out of sequence";
#endif

	  if ((Event.getType() == NONE)
	      || ((Event.getdt() > sorter->next_dt()) 
		  && (++_interactionRejectionCounter < rejectionLimit)))
	    {
	      this->fullUpdate(p1, p2);
	      return;
	    }

	  //Reset the rejection watchdog, we will run an interaction event now
	  _interactionRejectionCounter = 0;
		
#ifdef DYNAMO_DEBUG
	  if (Event.getdt() < 0)
	    derr << "Negative time " << Event.getdt() << std::endl;

	  if (boost::math::isnan(Event.getdt()))
	    M_throw() << "A NAN Interaction collision time has been found"
		      << Event.stringData(Sim);
	
	  if (Event.getdt() == HUGE_VAL)
	    M_throw() << "An infinite Interaction (not marked as NONE) collision time has been found\n"
		      << Event.stringData(Sim);
#endif
	
	  //Debug section
#ifdef DYNAMO_CollDebug
	  dout << "Interaction at t=" << (Event.getdt() + Sim->systemTime) / Sim->units.unitTime()
	       << "  ID1 " << ((p2.getID() < p2.getID()) ? p1.getID() : p2.getID())
	       << "  ID2 " << ((p2.getID() < p2.getID()) ? p2.getID() : p1.getID())
	       << "  dt " << Event.getdt()
	       << "  Type " << Event.getType() 
	       << std::endl;
#endif

	  Sim->systemTime += Event.getdt();
	
	  stream(Event.getdt());
	
	  //dynamics must be updated first
	  Sim->stream(Event.getdt());
	
	  Sim->interactions[Event.getInteractionID()]
	    ->runEvent(p1,p2,Event);

	  break;
	}
      case GLOBAL:
	{
	  //We don't stream the system for globals as neighbour lists
	  //optimise this (they dont need it).

	  //We also don't recheck Global events! (Check, some events might rely on this behavior)
	  Sim->globals[sorter->next_p2()]
	    ->runEvent(Sim->particles[sorter->next_ID()], sorter->next_dt());       	
	  break;	           
	}
      case LOCAL:
	{
	  Particle& part(Sim->particles[sorter->next_ID()]);

	  //Copy the FEL event
	  size_t localID = sorter->next_p2();

	  //Ready the next event in the FEL
	  sorter->popNextEvent();
	  sorter->update(sorter->next_ID());
	  sorter->sort();
	  lazyDeletionCleanup();

	  Sim->dynamics->updateParticle(part);
	  LocalEvent iEvent(Sim->locals[localID]->getEvent(part));

	  double next_dt = sorter->next_dt();

	  //Check the recalculated event is valid and not later than
	  //the next event in the queue
	  if ((iEvent.getType() == NONE)
	      || ((iEvent.getdt() > next_dt) && (++_localRejectionCounter < rejectionLimit)))
	    {
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
	
	  Sim->systemTime += iEvent.getdt();
	
	  stream(iEvent.getdt());
	
	  //dynamics must be updated first
	  Sim->stream(iEvent.getdt());
	
	  Sim->locals[localID]->runEvent(part, iEvent);	  
	  break;
	}
      case SYSTEM:
	{
	  Sim->systems[sorter->next_p2()]
	    ->runEvent();
	  //This saves the system events rebuilding themselves
	  rebuildSystemEvents();
	  break;
	}
      case RECALCULATE:
	{
	  //This is a special event type which requires that the
	  // events for this particle recalculated.
	  size_t ID = sorter->next_ID();
	  this->fullUpdate(Sim->particles[ID]);
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
    Particle& part1(Sim->particles[part.getID()]);
    Particle& part2(Sim->particles[id]);

    Sim->dynamics->updateParticle(part2);

    const IntEvent& eevent(Sim->getEvent(part1, part2));

    if (eevent.getType() != NONE)
      sorter->push(Event(eevent, eventCount[id]), part1.getID());
  }

  void 
  Scheduler::addLocalEvent(const Particle& part, 
			   const size_t& id) const
  {
    if (Sim->locals[id]->isInteraction(part))
      sorter->push(Sim->locals[id]->getEvent(part), part.getID());
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
