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

#include <dynamo/NparticleEventData.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/globals/PBCSentinel.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/simulation.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
GPBCSentinel::GPBCSentinel(dynamo::Simulation *nSim, const std::string &name)
    : Global(nSim, "PBCSentinel"), maxintdist(0) {
  globName = name;
  dout << "PBCSentinel Loaded" << std::endl;
}

GPBCSentinel::GPBCSentinel(const magnet::xml::Node &XML,
                           dynamo::Simulation *ptrSim)
    : Global(ptrSim, "PBCSentinel"), maxintdist(0) {
  GPBCSentinel::operator<<(XML);

  dout << "PBCSentinel Loaded" << std::endl;
}

void GPBCSentinel::initialise(size_t nID) {
  Global::initialise(nID);
  maxintdist = Sim->getLongestInteraction();
}

void GPBCSentinel::operator<<(const magnet::xml::Node &XML) {
  globName = XML.getAttribute("Name");
}

Event GPBCSentinel::getEvent(const Particle &part) const {
  return Event(part, Sim->dynamics->getPBCSentinelTime(part, maxintdist),
               GLOBAL, VIRTUAL, ID);
}

void GPBCSentinel::runEvent(Particle &part, const double dt) {
  Event iEvent(part, dt, GLOBAL, VIRTUAL, ID);

  Sim->systemTime += iEvent._dt;

  Sim->scheduler->stream(iEvent._dt);

  Sim->stream(iEvent._dt);

  NEventData EDat(ParticleEventData(part, *Sim->species[part], VIRTUAL));

  Sim->_sigParticleUpdate(EDat);

  for (shared_ptr<OutputPlugin> &Ptr : Sim->outputPlugins)
    Ptr->eventUpdate(iEvent, EDat);

  Sim->scheduler->fullUpdate(part);
}
} // namespace dynamo
