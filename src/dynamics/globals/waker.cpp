/*  DYNAMO:- Event driven molecular dynamics simulator 
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

#include "waker.hpp"
#include "globEvent.hpp"
#include "../NparticleEventData.hpp"
#include "../../base/is_simdata.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../../schedulers/scheduler.hpp"
#include "../units/units.hpp"

GWaker::GWaker(DYNAMO::SimData* nSim, const std::string& name, CRange* range, 
	       const double wt,const double wv):
  Global(range, nSim, "GWaker"),
  _wakeTime(wt),
  _wakeVelocity(wv)
{
  globName = name;
  I_cout() << "GWaker Loaded";
}

GWaker::GWaker(const XMLNode &XML, DYNAMO::SimData* ptrSim):
  Global(NULL, ptrSim, "GWaker")
{
  operator<<(XML);

  I_cout() << "GWaker Loaded";
}

void 
GWaker::initialise(size_t nID)
{
  ID=nID;
}

void 
GWaker::operator<<(const XMLNode& XML)
{
  range.set_ptr(CRange::loadClass(XML, Sim));

  try {
    globName = XML.getAttribute("Name");

    _wakeTime = Sim->dynamics.units().unitTime() * 
      boost::lexical_cast<double>(XML.getAttribute("WakeTime"));

    _wakeVelocity = Sim->dynamics.units().unitVelocity() * 
      boost::lexical_cast<double>(XML.getAttribute("WakeVelocity"));
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
  ParticleEventData EDat(part, Sim->dynamics.getSpecies(part), iEvent.getType());

  Vector newVel(Sim->normal_sampler(),Sim->normal_sampler(),Sim->normal_sampler());
  newVel *= _wakeVelocity / newVel.nrm();

  const_cast<Particle&>(part).getVelocity() = newVel;
  const_cast<Particle&>(part).setState(Particle::DYNAMIC);

  EDat.setDeltaKE(0.5 * EDat.getSpecies().getMass()
		  * (part.getVelocity().nrm2() 
		     - EDat.getOldVel().nrm2()));

  Sim->signalParticleUpdate(EDat);

  //Now we're past the event, update the scheduler and plugins
  Sim->ptrScheduler->fullUpdate(part);

  BOOST_FOREACH(magnet::ClonePtr<OutputPlugin> & Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(iEvent, EDat);

}

void 
GWaker::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Type") << "Waker"
      << xml::attr("Name") << globName
      << xml::attr("WakeVelocity") << _wakeVelocity / Sim->dynamics.units().unitVelocity()
      << xml::attr("WakeTime") << _wakeTime / Sim->dynamics.units().unitTime()
      << range;
}
