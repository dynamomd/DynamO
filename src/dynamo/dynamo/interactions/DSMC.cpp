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

#include <cmath>
#include <dynamo/2particleEventData.hpp>
#include <dynamo/BC/BC.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/interactions/DSMC.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <dynamo/ranges/IDRange.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/units/units.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
IDSMC::IDSMC(const magnet::xml::Node &XML, dynamo::Simulation *tmp)
    : ICapture(tmp, NULL) {
  IDSMC::operator<<(XML);
}

void IDSMC::initialise(size_t nID) {
  Interaction::initialise(nID);
  ICapture::initCaptureMap();
}

std::array<double, 4> IDSMC::getGlyphSize(size_t ID) const {
  return {{_length->getProperty(ID), 0, 0, 0}};
}

void IDSMC::operator<<(const magnet::xml::Node &XML) {
  Interaction::operator<<(XML);
  _length = Sim->_properties.getProperty(XML.getAttribute("Length"),
                                         Property::Units::Length());
  _e = Sim->_properties.getProperty(XML.getAttribute("Elasticity"),
                                    Property::Units::Dimensionless());
  ICapture::loadCaptureMap(XML);
}

double IDSMC::maxIntDist() const { return _length->getMaxValue(); }

Event IDSMC::getEvent(const Particle &p1, const Particle &p2) const {
#ifdef DYNAMO_DEBUG
  if (!Sim->dynamics->isUpToDate(p1))
    M_throw() << "Particle 1 is not up to date";

  if (!Sim->dynamics->isUpToDate(p2))
    M_throw() << "Particle 2 is not up to date";

  if (p1 == p2)
    M_throw() << "You shouldn't pass p1==p2 events to the interactions!";
#endif

  const double l = _length->getProperty(p1, p2);
  if (isCaptured(p1, p2)) {
    double dt = Sim->dynamics->SphereSphereOutRoot(p1, p2, l);
    return Event(p1, dt, INTERACTION, NBHOOD_OUT, ID, p2);
  } else {
    double dt = Sim->dynamics->SphereSphereInRoot(p1, p2, l);
    if (dt != std::numeric_limits<float>::infinity())
      return Event(p1, dt, INTERACTION, NBHOOD_IN, ID, p2);
  }

  return Event(p1, std::numeric_limits<float>::infinity(), INTERACTION, NONE,
               ID, p2);
}

PairEventData IDSMC::runEvent(Particle &p1, Particle &p2, Event iEvent) {
  PairEventData retval;

  ++Sim->eventCount;
  switch (iEvent._type) {
  case NBHOOD_IN:
    ICapture::add(p1, p2);
    return PairEventData(p1, p2, *Sim->species[p1], *Sim->species[p2], VIRTUAL);
  case NBHOOD_OUT:
    ICapture::remove(p1, p2);
    return PairEventData(p1, p2, *Sim->species[p1], *Sim->species[p2], VIRTUAL);
  case VIRTUAL:
    return PairEventData(p1, p2, *Sim->species[p1], *Sim->species[p2], VIRTUAL);
  default:
    M_throw() << "Unknown collision type";
  }
}

void IDSMC::outputXML(magnet::xml::XmlStream &XML) const {
  XML << magnet::xml::attr("Type") << "DSMC" << magnet::xml::attr("Length")
      << _length->getName() << magnet::xml::attr("Elasticity") << _e->getName()
      << magnet::xml::attr("Name") << intName << range;

  ICapture::outputCaptureMap(XML);
}

size_t IDSMC::captureTest(const Particle &p1, const Particle &p2) const {
  if (&(*(Sim->getInteraction(p1, p2))) != this)
    return false;

  const double l = _length->getProperty(p1, p2);
  return Sim->dynamics->sphereOverlap(p1, p2, l) > 0;
}

bool IDSMC::validateState(const Particle &p1, const Particle &p2,
                          bool textoutput) const {
  const double l = _length->getProperty(p1, p2);

  if (isCaptured(p1, p2)) {
    if (!Sim->dynamics->sphereOverlap(p1, p2, l)) {
      if (textoutput)
        derr << "Particle " << p1.getID() << " and Particle " << p2.getID()
             << " are registered as being closer than "
             << l / Sim->units.unitLength() << " but they are at a distance of "
             << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
             << std::endl;

      return true;
    }
  } else if (Sim->dynamics->sphereOverlap(p1, p2, l)) {
    if (textoutput)
      derr << "Particle " << p1.getID() << " and Particle " << p2.getID()
           << " are not registered as being closer than "
           << l / Sim->units.unitLength() << " but they are at a distance of "
           << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
           << std::endl;
    return true;
  }
  return false;
}
} // namespace dynamo
