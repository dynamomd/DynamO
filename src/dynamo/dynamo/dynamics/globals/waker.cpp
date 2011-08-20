/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
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

#include <dynamo/dynamics/globals/waker.hpp>
#include <dynamo/dynamics/globals/globEvent.hpp>
#include <dynamo/dynamics/NparticleEventData.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <dynamo/dynamics/liouvillean/liouvillean.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/dynamics/units/units.hpp>
#include <dynamo/dynamics/globals/neighbourList.hpp>
#include <magnet/xmlreader.hpp>

#ifdef DYNAMO_DEBUG 
#include <boost/math/special_functions/fpclassify.hpp>
#endif

namespace dynamo {
  GWaker::GWaker(dynamo::SimData* nSim, const std::string& name, CRange* range, 
		 const double wt,const double wv, std::string nblist):
    Global(nSim, "GWaker", range),
    _wakeTime(wt),
    _wakeVelocity(wv),
    _nblistName(nblist)
  {
    globName = name;
    dout << "GWaker Loaded" << std::endl;
  }

  GWaker::GWaker(const magnet::xml::Node& XML, dynamo::SimData* ptrSim):
    Global(ptrSim, "GWaker")
  {
    operator<<(XML);

    dout << "GWaker Loaded" << std::endl;
  }

  void 
  GWaker::initialise(size_t nID)
  {
    ID=nID;

    try {
      _NBListID = Sim->dynamics.getGlobal(_nblistName)->getID();
    }
    catch(std::exception& cxp)
      {
	M_throw() << "Failed while finding the neighbour list global.\n"
		  << "You must have a neighbour list for this waker event"
		  << cxp.what();
      }
  
    if (dynamic_cast<const GNeighbourList*>
	(Sim->dynamics.getGlobals()[_NBListID].get_ptr())
	== NULL)
      M_throw() << "The Global named SchedulerNBList is not a neighbour list!";

  }

  void 
  GWaker::operator<<(const magnet::xml::Node& XML)
  {
    range.set_ptr(CRange::getClass(XML, Sim));

    try {
      globName = XML.getAttribute("Name");

      _wakeTime = XML.getAttribute("WakeTime").as<double>() * Sim->dynamics.units().unitTime();

      _wakeVelocity = XML.getAttribute("WakeVelocity").as<double>() 
	* Sim->dynamics.units().unitVelocity();

      _nblistName = XML.getAttribute("NBList");
    }
    catch(...)
      {
	M_throw() << "Error loading GWaker";
      }
  }

  GlobalEvent
  GWaker::getEvent(const Particle& part) const
  {
    if (part.testState(Particle::DYNAMIC))
      return GlobalEvent(part,  HUGE_VAL, NONE, *this);
    else
      return GlobalEvent(part, _wakeTime, WAKEUP, *this);
  }

  void 
  GWaker::runEvent(const Particle& part, const double dt) const
  {
    GlobalEvent iEvent(getEvent(part));
    iEvent.setdt(dt); //We only trust the schedulers time, as we don't
    //track the motion of the system in Globals
  
#ifdef DYNAMO_DEBUG 
    if (boost::math::isnan(iEvent.getdt()))
      M_throw() << "A NAN Interaction collision time has been found"
		<< iEvent.stringData(Sim);
  
    if (iEvent.getdt() == HUGE_VAL)
      M_throw() << "An infinite Interaction (not marked as NONE) collision time has been found\n"
		<< iEvent.stringData(Sim);
#endif

    Sim->dSysTime += iEvent.getdt();
    
    Sim->ptrScheduler->stream(iEvent.getdt());
  
    Sim->dynamics.stream(iEvent.getdt());

    Sim->dynamics.getLiouvillean().updateParticle(part);

    //Here is where the particle goes to sleep or wakes
    ++Sim->eventCount;
  
    _neighbors = 0;
    //Grab a reference to the neighbour list
    const GNeighbourList& nblist(*static_cast<const GNeighbourList*>(Sim->dynamics.getGlobals()[_NBListID]
								     .get_ptr()));
    //Add the interaction events
    nblist.getParticleNeighbourhood(part, magnet::function::MakeDelegate(this, &GWaker::nblistCallback));  
  
    //  if (_neighbors < 10)
    //    {
    iEvent.addTime(Sim->freestreamAcc);      
    Sim->freestreamAcc = 0;

    ParticleEventData EDat(part, Sim->dynamics.getSpecies(part), iEvent.getType());
      
    Vector newVel(Sim->normal_sampler(),Sim->normal_sampler(),Sim->normal_sampler());
    newVel *= _wakeVelocity / newVel.nrm();
      
    const_cast<Particle&>(part).getVelocity() = newVel;
    const_cast<Particle&>(part).setState(Particle::DYNAMIC);
      
    EDat.setDeltaKE(0.5 * EDat.getSpecies().getMass(part.getID())
		    * (part.getVelocity().nrm2() 
		       - EDat.getOldVel().nrm2()));
      
    Sim->signalParticleUpdate(EDat);
      
    BOOST_FOREACH(magnet::ClonePtr<OutputPlugin> & Ptr, Sim->outputPlugins)
      Ptr->eventUpdate(iEvent, EDat);
    //    }
    //  else
    //    Sim->freestreamAcc += iEvent.getdt();

    //Now we're past the event, update the scheduler and plugins
    Sim->ptrScheduler->fullUpdate(part);
  }

  void 
  GWaker::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::tag("Global")
	<< magnet::xml::attr("Type") << "Waker"
	<< magnet::xml::attr("Name") << globName
	<< magnet::xml::attr("WakeVelocity") << _wakeVelocity / Sim->dynamics.units().unitVelocity()
	<< magnet::xml::attr("WakeTime") << _wakeTime / Sim->dynamics.units().unitTime()
	<< magnet::xml::attr("NBList") << _nblistName
	<< range
	<< magnet::xml::endtag("Global");
  }

  void 
  GWaker::nblistCallback(const Particle& part, const size_t& oid) const
  {
    Vector sep = part.getPosition() - Sim->particleList[oid].getPosition();
    Sim->dynamics.BCs().applyBC(sep);
    if (sep.nrm() < 2.01 * Sim->dynamics.units().unitLength())
      ++_neighbors;
  }
}
