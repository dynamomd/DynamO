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

#include <cmath> //for huge val
#include <dynamo/ranges/IDRangeNone.hpp>
#include <dynamo/schedulers/systemonly.hpp>
#include <dynamo/simulation.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
SSystemOnly::SSystemOnly(const magnet::xml::Node &XML,
                         dynamo::Simulation *const Sim)
    : Scheduler(Sim, "SystemOnlyScheduler", NULL) {
  dout << "System Events Only Scheduler Algorithm" << std::endl;
  operator<<(XML);
}

SSystemOnly::SSystemOnly(dynamo::Simulation *const Sim, FEL *ns)
    : Scheduler(Sim, "SystemOnlyScheduler", ns) {
  dout << "System Events Only Scheduler Algorithm" << std::endl;
}

void SSystemOnly::initialise() {
  dout << "Reinitialising on collision " << Sim->eventCount << std::endl;

  if (Sim->systems.empty())
    M_throw() << "A SystemOnlyScheduler used when there are no system events?";

  sorter->clear();
  sorter->init(Sim->N() + 1);
  rebuildSystemEvents();
}

void SSystemOnly::rebuildList() {
#ifdef DYNAMO_DEBUG
  initialise();
#else
  if (Sim->systems.empty())
    M_throw() << "A SystemOnlyScheduler used when there are no system events?";

  sorter->clear();
  sorter->init(Sim->N() + 1);
  rebuildSystemEvents();
#endif
}

void SSystemOnly::outputXML(magnet::xml::XmlStream &XML) const {
  XML << magnet::xml::attr("Type") << "SystemOnly" << magnet::xml::tag("Sorter")
      << *sorter << magnet::xml::endtag("Sorter");
}

std::unique_ptr<IDRange>
SSystemOnly::getParticleNeighbours(const Particle &) const {
  return std::unique_ptr<IDRange>(new IDRangeNone());
}

std::unique_ptr<IDRange>
SSystemOnly::getParticleNeighbours(const Vector &) const {
  return std::unique_ptr<IDRange>(new IDRangeNone());
}

std::unique_ptr<IDRange>
SSystemOnly::getParticleLocals(const Particle &) const {
  return std::unique_ptr<IDRange>(new IDRangeNone());
}
} // namespace dynamo
