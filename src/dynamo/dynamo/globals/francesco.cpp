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
#include <dynamo/globals/francesco.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/simulation.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
GFrancesco::GFrancesco(const magnet::xml::Node &XML, dynamo::Simulation *ptrSim)
    : Global(ptrSim, "Francesco") {
  operator<<(XML);
}

void GFrancesco::initialise(size_t nID) {
  Global::initialise(nID);
  _eventTimes.clear();
  _eventTimes.resize(Sim->particles.size(),
                     std::numeric_limits<float>::infinity());
  Sim->_sigParticleUpdate.connect<GFrancesco, &GFrancesco::particlesUpdated>(
      this);
}

double GFrancesco::generateTime() const {
  double dt;
  do
    dt = _dist(Sim->ranGenerator);
  while ((dt < 0) || (dt > (2 * _dist.mean())));
  return dt + Sim->systemTime;
}

void GFrancesco::particlesUpdated(const NEventData &PDat) {
  for (const PairEventData &pdat : PDat.L2partChanges) {
    _eventTimes[pdat.particle1_.getParticleID()] = generateTime();
    _eventTimes[pdat.particle2_.getParticleID()] = generateTime();
  }
}
void GFrancesco::outputXML(magnet::xml::XmlStream &XML) const {
  XML << magnet::xml::tag("Global") << magnet::xml::attr("Type") << "Francesco"
      << magnet::xml::attr("Name") << globName << magnet::xml::attr("MFT")
      << _dist.mean() / Sim->units.unitTime() << magnet::xml::attr("MFTstddev")
      << _dist.stddev() / Sim->units.unitTime() << magnet::xml::attr("Velocity")
      << _vel / Sim->units.unitVelocity() << range
      << magnet::xml::endtag("Global");
}

void GFrancesco::operator<<(const magnet::xml::Node &XML) {
  globName = XML.getAttribute("Name");
  const double MFT =
      XML.getAttribute("MFT").as<double>() * Sim->units.unitTime();
  const double MFTstddev =
      XML.getAttribute("MFTstddev").as<double>() * Sim->units.unitTime();
  _dist = std::normal_distribution<>(MFT, MFTstddev);
  _vel = XML.getAttribute("Velocity").as<double>() * Sim->units.unitVelocity();
  range = shared_ptr<IDRange>(IDRange::getClass(XML.getNode("IDRange"), Sim));
}

Event GFrancesco::getEvent(const Particle &part) const {
  return Event(part, _eventTimes[part.getID()] - Sim->systemTime, GLOBAL,
               GAUSSIAN, ID);
}

void GFrancesco::runEvent(Particle &part, const double) {
  const double dt = _eventTimes[part.getID()] - Sim->systemTime;
  _eventTimes[part.getID()] = std::numeric_limits<float>::infinity();
  Event iEvent(part, dt, GLOBAL, GAUSSIAN, ID);

  Sim->systemTime += dt;
  Sim->ptrScheduler->stream(dt);
  Sim->stream(dt);

  Sim->dynamics->updateParticle(part);
  NEventData EDat(ParticleEventData(part, *Sim->species[part], GAUSSIAN));

  // Kill the rotational motion
  Sim->dynamics->getRotData(part).angularVelocity = Vector{0, 0, 0};
  // Reassign the linear motion
  part.getVelocity() = _vel * (Sim->dynamics->getRotData(part).orientation *
                               magnet::math::Quaternion::initialDirector());

  Sim->_sigParticleUpdate(EDat);
  for (shared_ptr<OutputPlugin> &Ptr : Sim->outputPlugins)
    Ptr->eventUpdate(iEvent, EDat);
  Sim->ptrScheduler->fullUpdate(part);
}
} // namespace dynamo
