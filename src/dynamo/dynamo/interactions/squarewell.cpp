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

#include <dynamo/BC/BC.hpp>
#include <dynamo/interactions/squarewell.hpp>

#include <cmath>
#include <dynamo/2particleEventData.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/globals/global.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <dynamo/particle.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/species/species.hpp>
#include <dynamo/units/units.hpp>
#include <iomanip>
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
ISquareWell::ISquareWell(const magnet::xml::Node &XML, dynamo::Simulation *tmp)
    : ICapture(tmp, NULL) // A temporary value!
{
  ISquareWell::operator<<(XML);
}

void ISquareWell::operator<<(const magnet::xml::Node &XML) {
  Interaction::operator<<(XML);
  _diameter = Sim->_properties.getProperty(XML.getAttribute("Diameter"),
                                           Property::Units::Length());
  _lambda = Sim->_properties.getProperty(XML.getAttribute("Lambda"),
                                         Property::Units::Dimensionless());
  _wellDepth = Sim->_properties.getProperty(XML.getAttribute("WellDepth"),
                                            Property::Units::Energy());

  if (XML.hasAttribute("Elasticity"))
    _e = Sim->_properties.getProperty(XML.getAttribute("Elasticity"),
                                      Property::Units::Dimensionless());
  else
    _e = Sim->_properties.getProperty(1.0, Property::Units::Dimensionless());

  ICapture::loadCaptureMap(XML);
}

std::array<double, 4> ISquareWell::getGlyphSize(size_t ID) const {
  return {{_diameter->getProperty(ID), 0, 0, 0}};
}

double ISquareWell::getExcludedVolume(size_t ID) const {
  double diam = _diameter->getProperty(ID);
  return diam * diam * diam * M_PI / 6.0;
}

double ISquareWell::maxIntDist() const {
  return _diameter->getMaxValue() * _lambda->getMaxValue();
}

void ISquareWell::initialise(size_t nID) {
  Interaction::initialise(nID);
  ICapture::initCaptureMap();
}

size_t ISquareWell::captureTest(const Particle &p1, const Particle &p2) const {
  if (&(*(Sim->getInteraction(p1, p2))) != this)
    return false;

  const double d = _diameter->getProperty(p1, p2);
  const double l = _lambda->getProperty(p1, p2);

#ifdef DYNAMO_DEBUG
  if (Sim->dynamics->sphereOverlap(p1, p2, d))
    derr << "Warning! Two particles might be overlapping"
         << "Overlap is "
         << Sim->dynamics->sphereOverlap(p1, p2, d) / Sim->units.unitLength()
         << "\nd = " << d / Sim->units.unitLength() << std::endl;
#endif

  return Sim->dynamics->sphereOverlap(p1, p2, l * d) > 0;
}

Event ISquareWell::getEvent(const Particle &p1, const Particle &p2) const {
#ifdef DYNAMO_DEBUG
  if (!Sim->dynamics->isUpToDate(p1))
    M_throw() << "Particle 1 is not up to date";

  if (!Sim->dynamics->isUpToDate(p2))
    M_throw() << "Particle 2 is not up to date";

  if (p1 == p2)
    M_throw() << "You shouldn't pass p1==p2 events to the interactions!";
#endif

  const double d = _diameter->getProperty(p1, p2);
  const double l = _lambda->getProperty(p1, p2);

  Event retval(p1, std::numeric_limits<float>::infinity(), INTERACTION, NONE,
               ID, p2);

  if (isCaptured(p1, p2)) {
    double dt = Sim->dynamics->SphereSphereInRoot(p1, p2, d);
    if (dt != std::numeric_limits<float>::infinity())
      retval = Event(p1, dt, INTERACTION, CORE, ID, p2);

    dt = Sim->dynamics->SphereSphereOutRoot(p1, p2, l * d);
    if (retval._dt > dt)
      retval = Event(p1, dt, INTERACTION, STEP_OUT, ID, p2);
  } else {
    double dt = Sim->dynamics->SphereSphereInRoot(p1, p2, l * d);

    if (dt != std::numeric_limits<float>::infinity())
      retval = Event(p1, dt, INTERACTION, STEP_IN, ID, p2);
  }

  return retval;
}

PairEventData ISquareWell::runEvent(Particle &p1, Particle &p2, Event iEvent) {
  ++Sim->eventCount;

  const double d = _diameter->getProperty(p1, p2);
  const double d2 = d * d;
  const double e = _e->getProperty(p1, p2);
  const double l = _lambda->getProperty(p1, p2);
  const double ld2 = d * l * d * l;
  const double wd = _wellDepth->getProperty(p1, p2);

  PairEventData retVal;
  switch (iEvent._type) {
  case CORE:
    return Sim->dynamics->SmoothSpheresColl(iEvent, e, d2, CORE);
  case STEP_IN: {
    PairEventData retVal = Sim->dynamics->SphereWellEvent(iEvent, wd, ld2, 1);
    if (retVal.getType() != BOUNCE)
      ICapture::add(p1, p2);
    return retVal;
  }
  case STEP_OUT: {
    PairEventData retVal = Sim->dynamics->SphereWellEvent(iEvent, -wd, ld2, 0);
    if (retVal.getType() != BOUNCE)
      ICapture::remove(p1, p2);
    return retVal;
  }
  default:
    M_throw() << "Unknown collision type";
  }
}

bool ISquareWell::validateState(const Particle &p1, const Particle &p2,
                                bool textoutput) const {
  const double d = _diameter->getProperty(p1, p2);
  const double l = _lambda->getProperty(p1, p2);

  if (isCaptured(p1, p2)) {
    if (!Sim->dynamics->sphereOverlap(p1, p2, l * d)) {
      if (textoutput)
        derr << "Particle " << p1.getID() << " and Particle " << p2.getID()
             << " registered as being inside the well at "
             << l * d / Sim->units.unitLength()
             << " but they are at a distance of "
             << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
             << std::endl;

      return true;
    }

    if (Sim->dynamics->sphereOverlap(p1, p2, d)) {
      if (textoutput)
        derr << "Particle " << p1.getID() << " and Particle " << p2.getID()
             << " are inside the well with an inner hard core at "
             << d / Sim->units.unitLength() << " but they are at a distance of "
             << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
             << std::endl;

      return true;
    }
  } else if (Sim->dynamics->sphereOverlap(p1, p2, l * d)) {
    if (textoutput)
      derr << "Particle " << p1.getID() << " and Particle " << p2.getID()
           << " are registered as being outside the well at a distance of "
           << l * d / Sim->units.unitLength()
           << " but they are at a distance of "
           << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
           << std::endl;

    return true;
  }

  return false;
}

void ISquareWell::outputXML(magnet::xml::XmlStream &XML) const {
  XML << magnet::xml::attr("Type") << "SquareWell"
      << magnet::xml::attr("Diameter") << _diameter->getName()
      << magnet::xml::attr("Elasticity") << _e->getName()
      << magnet::xml::attr("Lambda") << _lambda->getName()
      << magnet::xml::attr("WellDepth") << _wellDepth->getName()
      << magnet::xml::attr("Name") << intName << *range;

  ICapture::outputCaptureMap(XML);
}

double ISquareWell::getInternalEnergy() const {
  // Once the capture maps are loaded just iterate through that determining
  // energies
  double Energy = 0.0;
  for (const ICapture::value_type &IDs : *this)
    Energy += getInternalEnergy(Sim->particles[IDs.first.first],
                                Sim->particles[IDs.first.second]);
  return Energy;
}

double ISquareWell::getInternalEnergy(const Particle &p1,
                                      const Particle &p2) const {
  return -_wellDepth->getProperty(p1, p2) * isCaptured(p1, p2);
}
} // namespace dynamo
