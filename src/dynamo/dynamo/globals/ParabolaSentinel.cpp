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

#include <dynamo/globals/ParabolaSentinel.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {

  GParabolaSentinel::GParabolaSentinel(dynamo::Simulation* nSim, const std::string& name):
    Global(nSim, "ParabolaSentinel")
  {
    globName = name;
    dout << "ParabolaSentinel Loaded" << std::endl;
  }

  Event
  GParabolaSentinel::getEvent(const Particle& part) const
  {
    Sim->dynamics->updateParticle(Sim->particles[part.getID()]);

    return Event(part, Sim->dynamics->getParabolaSentinelTime(part), GLOBAL, RECALCULATE_PARABOLA, ID);
  }

  void 
  GParabolaSentinel::runEvent(Particle& part, const double)
  {
    Sim->dynamics->updateParticle(part);

    Event iEvent = getEvent(part);

    iEvent._type = VIRTUAL;

    if (iEvent._dt == std::numeric_limits<float>::infinity())
      {
	//We've numerically drifted slightly passed the parabola, so
	//just reschedule the particles events, no need to enforce anything
	Sim->ptrScheduler->fullUpdate(part);
	return;
      }

#ifdef DYNAMO_DEBUG 
    if (std::isnan(iEvent._dt))
      M_throw() << "A NAN Interaction collision time has been found when recalculating this global";
#endif

    Sim->systemTime += iEvent._dt;
    
    Sim->ptrScheduler->stream(iEvent._dt);
  
    Sim->stream(iEvent._dt);

    NEventData EDat = Sim->dynamics->enforceParabola(part);
  
    Sim->_sigParticleUpdate(EDat);

    for (shared_ptr<OutputPlugin> & Ptr : Sim->outputPlugins)
      Ptr->eventUpdate(iEvent, EDat);

    Sim->ptrScheduler->fullUpdate(part);
  }
}
