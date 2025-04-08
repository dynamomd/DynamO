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
#include <dynamo/outputplugins/outputplugin.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/systems/tHalt.hpp>
#include <dynamo/units/units.hpp>

namespace dynamo {
SystHalt::SystHalt(dynamo::Simulation *nSim, double ndt, std::string nName)
    : System(nSim) {
  dt = ndt * Sim->units.unitTime();

  sysName = nName;

  dout << "System halt set for " << ndt << std::endl;

  type = VIRTUAL;
}

NEventData SystHalt::runEvent() {
  Sim->nextPrintEvent = Sim->endEventCount = Sim->eventCount;
  return NEventData();
}

void SystHalt::initialise(size_t nID) { ID = nID; }

void SystHalt::setdt(double ndt) { dt = ndt * Sim->units.unitTime(); }

void SystHalt::increasedt(double ndt) { dt += ndt * Sim->units.unitTime(); }
} // namespace dynamo
