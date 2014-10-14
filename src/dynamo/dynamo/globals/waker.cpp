/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>
    Copyright (C) 2011  Sebastian Gonzalez <tsuresuregusa@gmail.com>

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

#include <dynamo/globals/waker.hpp>
#include <dynamo/BC/BC.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/globals/neighbourList.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
  GWaker::GWaker(dynamo::Simulation* nSim, const std::string& name, IDRange* range, 
		 const double wt,const double wv, std::string nblist):
    Global(nSim, "GWaker", range),
    _wakeTime(wt),
    _wakeVelocity(wv),
    _nblistName(nblist)
  {
    globName = name;
    dout << "GWaker Loaded" << std::endl;
  }

  GWaker::GWaker(const magnet::xml::Node& XML, dynamo::Simulation* ptrSim):
    Global(ptrSim, "GWaker")
  {
    operator<<(XML);

    dout << "GWaker Loaded" << std::endl;
  }

  void 
  GWaker::initialise(size_t nID)
  {
    Global::initialise(nID);

    try {
      _NBListID = Sim->globals[_nblistName]->getID();
    }
    catch(std::exception& cxp)
      {
	M_throw() << "Failed while finding the neighbour list global.\n"
		  << "You must have a neighbour list for this waker event"
		  << cxp.what();
      }
  
    if (!std::dynamic_pointer_cast<GNeighbourList>(Sim->globals[_NBListID]))
      M_throw() << "The Global named SchedulerNBList is not a neighbour list!";

  }

  void 
  GWaker::operator<<(const magnet::xml::Node& XML)
  {
    range = shared_ptr<IDRange>(IDRange::getClass(XML, Sim));

    try {
      globName = XML.getAttribute("Name");

      _wakeTime = XML.getAttribute("WakeTime").as<double>() * Sim->units.unitTime();

      _wakeVelocity = XML.getAttribute("WakeVelocity").as<double>() 
	* Sim->units.unitVelocity();

      _nblistName = XML.getAttribute("NBList");
    }
    catch(...)
      {
	M_throw() << "Error loading GWaker";
      }
  }

  Event
  GWaker::getEvent(const Particle& part) const
  {
    if (part.testState(Particle::DYNAMIC))
      return Event(part, HUGE_VAL, GLOBAL, NONE, ID);
    else
      return Event(part, _wakeTime, GLOBAL, WAKEUP, ID);
  }

  void 
  GWaker::runEvent(Particle& part, const double dt)
  {
    Event iEvent = getEvent(part);
    iEvent._dt = dt; //We only trust the schedulers time, as we don't
    //track the motion of the system in Globals
  
#ifdef DYNAMO_DEBUG 
    if (std::isnan(iEvent._dt))
      M_throw() << "A NAN Interaction collision time has been found"
		<< iEvent;
  
    if (iEvent._dt == HUGE_VAL)
      M_throw() << "An infinite Interaction (not marked as NONE) collision time has been found\n"
		<< iEvent;
#endif

    Sim->systemTime += iEvent._dt;
    
    Sim->ptrScheduler->stream(iEvent._dt);
  
    Sim->stream(iEvent._dt);

    Sim->dynamics->updateParticle(part);

    //Here is where the particle goes to sleep or wakes
    ++Sim->eventCount;
  
    _neighbors = 0;

    //Add the interaction events
    std::unique_ptr<IDRange> ids(Sim->ptrScheduler->getParticleNeighbours(part));
    for (const size_t& id1 : *ids)
      nblistCallback(part, id1);
  
    ParticleEventData EDat(part, *Sim->species(part), iEvent._type);
    
    std::normal_distribution<> norm_dist;
    Vector newVel{norm_dist(Sim->ranGenerator), norm_dist(Sim->ranGenerator), norm_dist(Sim->ranGenerator)};
    newVel *= _wakeVelocity / newVel.nrm();
      
    part.getVelocity() = newVel;
    part.setState(Particle::DYNAMIC);
      
    Sim->_sigParticleUpdate(EDat);
      
    for (shared_ptr<OutputPlugin> & Ptr : Sim->outputPlugins)
      Ptr->eventUpdate(iEvent, EDat);

    //Now we're past the event, update the scheduler and plugins
    Sim->ptrScheduler->fullUpdate(part);
  }

  void 
  GWaker::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::tag("Global")
	<< magnet::xml::attr("Type") << "Waker"
	<< magnet::xml::attr("Name") << globName
	<< magnet::xml::attr("WakeVelocity") << _wakeVelocity / Sim->units.unitVelocity()
	<< magnet::xml::attr("WakeTime") << _wakeTime / Sim->units.unitTime()
	<< magnet::xml::attr("NBList") << _nblistName
	<< range
	<< magnet::xml::endtag("Global");
  }

  void 
  GWaker::nblistCallback(const Particle& part, const size_t& oid) const
  {
    Vector sep = part.getPosition() - Sim->particles[oid].getPosition();
    Sim->BCs->applyBC(sep);
    if (sep.nrm() < 2.01 * Sim->units.unitLength())
      ++_neighbors;
  }
}
