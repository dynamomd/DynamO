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
#include <dynamo/globals/globEvent.hpp>
#include <dynamo/BC/BC.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/simdata.hpp>
#include <dynamo/liouvillean/liouvillean.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/globals/neighbourList.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <magnet/xmlreader.hpp>

#ifdef DYNAMO_DEBUG 
#include <boost/math/special_functions/fpclassify.hpp>
#endif

namespace dynamo {
  GWaker::GWaker(dynamo::SimData* nSim, const std::string& name, Range* range, 
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
      _NBListID = Sim->globals[_nblistName].getID();
    }
    catch(std::exception& cxp)
      {
	M_throw() << "Failed while finding the neighbour list global.\n"
		  << "You must have a neighbour list for this waker event"
		  << cxp.what();
      }
  
    if (!std::tr1::dynamic_pointer_cast<GNeighbourList>(Sim->globals[_NBListID]))
      M_throw() << "The Global named SchedulerNBList is not a neighbour list!";

  }

  void 
  GWaker::operator<<(const magnet::xml::Node& XML)
  {
    range = shared_ptr<Range>(Range::getClass(XML, Sim));

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

  GlobalEvent
  GWaker::getEvent(const Particle& part) const
  {
    if (part.testState(Particle::DYNAMIC))
      return GlobalEvent(part,  HUGE_VAL, NONE, *this);
    else
      return GlobalEvent(part, _wakeTime, WAKEUP, *this);
  }

  void 
  GWaker::runEvent(Particle& part, const double dt) const
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
  
    Sim->stream(iEvent.getdt());

    Sim->liouvillean->updateParticle(part);

    //Here is where the particle goes to sleep or wakes
    ++Sim->eventCount;
  
    _neighbors = 0;

    //Add the interaction events
    Sim->ptrScheduler->getParticleNeighbourhood(part, magnet::function::MakeDelegate(this, &GWaker::nblistCallback));
  
    //  if (_neighbors < 10)
    //    {
    iEvent.addTime(Sim->freestreamAcc);      
    Sim->freestreamAcc = 0;

    ParticleEventData EDat(part, Sim->species[part], iEvent.getType());
      
    Vector newVel(Sim->normal_sampler(),Sim->normal_sampler(),Sim->normal_sampler());
    newVel *= _wakeVelocity / newVel.nrm();
      
    part.getVelocity() = newVel;
    part.setState(Particle::DYNAMIC);
      
    EDat.setDeltaKE(0.5 * EDat.getSpecies().getMass(part.getID())
		    * (part.getVelocity().nrm2() 
		       - EDat.getOldVel().nrm2()));
      
    Sim->signalParticleUpdate(EDat);
      
    BOOST_FOREACH(shared_ptr<OutputPlugin> & Ptr, Sim->outputPlugins)
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
	<< magnet::xml::attr("WakeVelocity") << _wakeVelocity / Sim->units.unitVelocity()
	<< magnet::xml::attr("WakeTime") << _wakeTime / Sim->units.unitTime()
	<< magnet::xml::attr("NBList") << _nblistName
	<< range
	<< magnet::xml::endtag("Global");
  }

  void 
  GWaker::nblistCallback(const Particle& part, const size_t& oid) const
  {
    Vector sep = part.getPosition() - Sim->particleList[oid].getPosition();
    Sim->BCs->applyBC(sep);
    if (sep.nrm() < 2.01 * Sim->units.unitLength())
      ++_neighbors;
  }
}
