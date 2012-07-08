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

#include <dynamo/globals/PBCSentinel.hpp>
#include <dynamo/globals/globEvent.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <magnet/xmlreader.hpp>

#ifdef DYNAMO_DEBUG 
#include <boost/math/special_functions/fpclassify.hpp>
#endif

namespace dynamo {
  GPBCSentinel::GPBCSentinel(dynamo::Simulation* nSim, const std::string& name):
    Global(nSim, "PBCSentinel"),
    maxintdist(0)
  {
    globName = name;
    dout << "PBCSentinel Loaded" << std::endl;
  }

  GPBCSentinel::GPBCSentinel(const magnet::xml::Node& XML, dynamo::Simulation* ptrSim):
    Global(ptrSim, "PBCSentinel"),
    maxintdist(0)
  {
    operator<<(XML);

    dout << "PBCSentinel Loaded" << std::endl;
  }

  void 
  GPBCSentinel::initialise(size_t nID)
  {
    ID=nID;
  
    maxintdist = Sim->getLongestInteraction();
  }

  void 
  GPBCSentinel::operator<<(const magnet::xml::Node& XML)
  {
    globName = XML.getAttribute("Name");	
  }

  GlobalEvent 
  GPBCSentinel::getEvent(const Particle& part) const
  {
    return GlobalEvent(part, Sim->dynamics
		       ->getPBCSentinelTime(part, maxintdist),
		       VIRTUAL, *this);
  }

  void 
  GPBCSentinel::runEvent(Particle& part, const double dt) const
  {
    GlobalEvent iEvent(part, dt, VIRTUAL, *this);

    Sim->systemTime += iEvent.getdt();
    
    Sim->ptrScheduler->stream(iEvent.getdt());
  
    Sim->stream(iEvent.getdt());

    NEventData EDat(ParticleEventData(part, *Sim->species[part], VIRTUAL));

    Sim->signalParticleUpdate(EDat);

    BOOST_FOREACH(shared_ptr<OutputPlugin> & Ptr, Sim->outputPlugins)
      Ptr->eventUpdate(iEvent, EDat);

    Sim->ptrScheduler->fullUpdate(part);
  }
}
