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

#include <dynamo/2particleEventData.hpp>
#include <dynamo/BC/BC.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/dynamics/compression.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/interactions/dumbbells.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <dynamo/ranges/IDRange.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/species/sphericalTop.hpp>
#include <dynamo/units/units.hpp>
#include <magnet/intersection/offcentre_spheres.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
IDumbbells::IDumbbells(const magnet::xml::Node &XML, dynamo::Simulation *tmp)
    : ICapture(tmp, NULL),
      _unusedDimension(std::numeric_limits<size_t>::max()) {
  operator<<(XML);
}

void IDumbbells::initialise(size_t nID) {
  Interaction::initialise(nID);
  ICapture::initCaptureMap();
}

std::array<double, 4> IDumbbells::getGlyphSize(size_t ID) const {

  return {{
      _compositeData.front()._diam->getProperty(ID), // DiamA
      _compositeData.back()._diam->getProperty(ID),  // DiamB
      +_compositeData.front()._offset.nrm(),         // LA
      -_compositeData.back()._offset.nrm()           // LB
  }};
}

double IDumbbells::getExcludedVolume(size_t ID) const {
  double vol = 0;
  for (const auto &p : _compositeData) {
    const double diam = p._diam->getProperty(ID);
    vol += diam * diam * diam * M_PI / 6.0;
  }
  return vol;
}

void IDumbbells::operator<<(const magnet::xml::Node &XML) {
  Interaction::operator<<(XML);
  _e = Sim->_properties.getProperty(XML.getAttribute("Elasticity"),
                                    Property::Units::Dimensionless());

  if (XML.hasAttribute("UnusedDimension"))
    _unusedDimension = XML.getAttribute("UnusedDimension").as<size_t>();

  for (magnet::xml::Node node = XML.findNode("Sphere"); node.valid(); ++node) {
    Vector offset;
    offset << node.getNode("Offset");
    _compositeData.push_back(Sphere(
        offset * Sim->units.unitLength(),
        Sim->_properties.getProperty(node.getAttribute("Diameter"),
                                     Property::Units::Length()),
        Sim->_properties.getProperty(0, Property::Units::Energy()),
        Sim->_properties.getProperty(0, Property::Units::Dimensionless())));
  }

  if (_compositeData.empty())
    M_throw() << "Interaction \"" << intName
              << "\" missing Sphere tags (requires at least 1)\nXML Path:"
              << XML.getPath();

  ICapture::loadCaptureMap(XML);
}

double IDumbbells::maxIntDist() const {
  double maxdist = 0;
  for (const auto &p : _compositeData) {
    const double diam = p._diam->getMaxValue();
    const double L = p._offset.nrm();
    maxdist = std::max(diam + 2 * L, maxdist);
  }
  return maxdist;
}

double IDumbbells::maxIntDist(size_t p1, size_t p2) const {
  double maxdist1 = 0, maxdist2 = 0;
  for (const auto &p : _compositeData) {
    const double diam1 = p._diam->getProperty(p1);
    const double L = p._offset.nrm();
    const double diam2 = p._diam->getProperty(p2);
    maxdist1 = std::max(0.5 * diam1 + L, maxdist1);
    maxdist2 = std::max(0.5 * diam2 + L, maxdist2);
  }
  return maxdist1 + maxdist2;
}

Event IDumbbells::getEvent(const Particle &p1, const Particle &p2) const {
#ifdef DYNAMO_DEBUG
  if (!Sim->dynamics->isUpToDate(p1))
    M_throw() << "Particle 1 is not up to date";

  if (!Sim->dynamics->isUpToDate(p2))
    M_throw() << "Particle 2 is not up to date";

  if (p1 == p2)
    M_throw() << "You shouldn't pass p1==p2 events to the interactions!";
#endif

  const double max_dist = maxIntDist();

  if (!isCaptured(p1, p2)) {
    double dt = Sim->dynamics->SphereSphereInRoot(p1, p2, max_dist);
    return Event(p1, dt, INTERACTION,
                 (dt != std::numeric_limits<float>::infinity()) ? NBHOOD_IN
                                                                : NONE,
                 ID, p2);
  }

  Vector r12 = p1.getPosition() - p2.getPosition();
  Vector v12 = p1.getVelocity() - p2.getVelocity();
  Sim->BCs->applyBC(r12, v12);
  const auto &angv1 = Sim->dynamics->getRotData(p1).angularVelocity;
  const auto &angv2 = Sim->dynamics->getRotData(p2).angularVelocity;

  double growthrate = 0;
  if (std::dynamic_pointer_cast<DynCompression>(Sim->dynamics))
    growthrate = std::static_pointer_cast<DynCompression>(Sim->dynamics)
                     ->getGrowthRate();

  // Run this to determine when the spheres no longer intersect
  const double t_max = Sim->dynamics->SphereSphereOutRoot(p1, p2, max_dist);

  std::pair<bool, double> current(false,
                                  std::numeric_limits<float>::infinity());
  // If the bounding spheres never stop intersecting, we need to
  // establish an upper time to search for events as the intersection
  // routine needs a maximum upper bound for the distance. If we
  // reach this time, we recalculate for events from then.
  if (t_max == std::numeric_limits<float>::infinity())
    current.second = 1.0;

  for (auto it1 = _compositeData.begin(); it1 != _compositeData.end(); ++it1)
    for (auto it2 = _compositeData.begin(); it2 != _compositeData.end();
         ++it2) {
      const double t_max_current = std::min(t_max, current.second);
      const Vector u1 =
          Sim->dynamics->getRotData(p1).orientation * it1->_offset;
      const Vector u2 =
          Sim->dynamics->getRotData(p2).orientation * it2->_offset;
      const double Diam1 = it1->_diam->getProperty(p1);
      const double Diam2 = it2->_diam->getProperty(p2);

      // const double max_growth_factor = 1 + std::max(growthrate *
      // Sim->systemTime, growthrate * (Sim->systemTime + t_max_current)); dout
      // << "max_dist=" << max_dist << std::endl; dout << "max_growth_factor="
      // << max_growth_factor <<std::endl; dout << "r12=" << r12.toString()
      // <<std::endl; dout << "v12=" << v12.toString() <<std::endl; dout <<
      // "angv1=" << angv1.toString() <<std::endl; dout << "angv2=" <<
      // angv2.toString() <<std::endl; dout << "u1=" << Vector(director1 *
      // L1).toString() << std::endl; dout << "u2=" << Vector(director2 *
      // L2).toString() << std::endl; dout << "Diam1=" << Diam1 << std::endl;
      // dout << "Diam2=" << Diam1 << std::endl;
      // dout << "max_dist_eff=" << max_dist * max_growth_factor << std::endl;
      // dout << "t=" << Sim->systemTime << std::endl;
      // dout << "growthrate=" << growthrate << std::endl;
      // dout << "tmin=" << 0 << std::endl;
      // dout << "tmax=" << t_max_current << std::endl;
      magnet::intersection::detail::OffcentreSpheresOverlapFunction f(
          r12, v12, angv1, angv2, u1, u2, Diam1, Diam2, max_dist,
          Sim->systemTime, growthrate, 0, t_max_current);

      std::pair<bool, double> test = f.nextEvent();
      if (test.second < current.second)
        current = test;
    }

  // Check if they miss each other
  if (current.second == std::numeric_limits<float>::infinity())
    return Event(p1, t_max, INTERACTION, NBHOOD_OUT, ID, p2);

  // Something happens in the time interval
  return Event(p1, current.second, INTERACTION, current.first ? CORE : VIRTUAL,
               ID, p2);
}

PairEventData IDumbbells::runEvent(Particle &p1, Particle &p2, Event iEvent) {
  switch (iEvent._type) {
  case CORE: {
    shared_ptr<SpSphericalTop> sp1 =
        std::dynamic_pointer_cast<SpSphericalTop>(Sim->species[p1]);
    shared_ptr<SpSphericalTop> sp2 =
        std::dynamic_pointer_cast<SpSphericalTop>(Sim->species[p2]);

    if (!sp1 || !sp2)
      M_throw() << "Could not find the intertia of one of the particles "
                   "undergoing an interaction";

    const Vector angvel1 = Sim->dynamics->getRotData(p1).angularVelocity;
    const Vector angvel2 = Sim->dynamics->getRotData(p2).angularVelocity;
    const double m1 = sp1->getMass(p1.getID());
    const double m2 = sp2->getMass(p2.getID());
    const double I1 = sp1->getScalarMomentOfInertia(p1.getID());
    const double I2 = sp2->getScalarMomentOfInertia(p2.getID());

    const double max_dist = maxIntDist(p1, p2);

    PairEventData retval = PairEventData(p1, p2, *sp1, *sp2, CORE);
    Sim->BCs->applyBC(retval.rij, retval.vijold);

    double growthrate = 0;
    if (std::dynamic_pointer_cast<DynCompression>(Sim->dynamics))
      growthrate = std::static_pointer_cast<DynCompression>(Sim->dynamics)
                       ->getGrowthRate();

    const double growthfactor = 1 + growthrate * Sim->systemTime;

    double fcurrent = std::numeric_limits<float>::infinity();
    double d1, d2;
    Vector l1, l2;
    d1 = d2 = 0;
    for (auto it1 = _compositeData.begin(); it1 != _compositeData.end(); ++it1)
      for (auto it2 = _compositeData.begin(); it2 != _compositeData.end();
           ++it2) {
        const Vector u1 =
            Sim->dynamics->getRotData(p1).orientation * it1->_offset;
        const Vector u2 =
            Sim->dynamics->getRotData(p2).orientation * it2->_offset;
        const double Diam1 = it1->_diam->getProperty(p1);
        const double Diam2 = it2->_diam->getProperty(p2);

        magnet::intersection::detail::OffcentreSpheresOverlapFunction f(
            retval.rij, retval.vijold, angvel1, angvel2, u1, u2, Diam1, Diam2,
            max_dist, Sim->systemTime, growthrate, 0, 0);
        std::array<double, 2> fval = f.eval<0, 2>();
        if ((fval[1] < 0) && (fval[0] < fcurrent)) {
          fcurrent = fval[0];
          d1 = Diam1;
          d2 = Diam2;
          l1 = u1;
          l2 = u2;
        }
      }

    // If no particles satisfy the collision condition, its a
    // numerical error so just return a virtual event
    if (fcurrent == std::numeric_limits<float>::infinity())
      return PairEventData(p1, p2, *Sim->species[p1], *Sim->species[p2],
                           VIRTUAL);

    ++Sim->eventCount;

    const Vector u1 = l1 * growthfactor;
    const Vector u2 = l2 * growthfactor;
    Vector nhat = retval.rij + u1 - u2;
    nhat /= nhat.nrm();
    const Vector r1 = u1 - nhat * 0.5 * d1 * growthfactor;
    const Vector r2 = u2 + nhat * 0.5 * d2 * growthfactor;
    const Vector vc12 = retval.vijold + (angvel1 ^ r1) - (angvel2 ^ r2) +
                        growthrate * (u1 - u2 - nhat * (d1 + d2) * 0.5);
    const double e = _e->getProperty(p1, p2);
    const double J = (1 + e) * (nhat | vc12) /
                     ((1 / m1) + (1 / m2) +
                      (nhat | ((1 / I1) * ((u1 ^ nhat) ^ u1) +
                               (1 / I2) * ((u2 ^ nhat) ^ u2))));

    retval.rvdot = (retval.rij | retval.vijold);
    retval.impulse = J * nhat;

    p1.getVelocity() -= retval.impulse / m1;
    p2.getVelocity() += retval.impulse / m2;
    Sim->dynamics->getRotData(p1).angularVelocity -= (r1 ^ retval.impulse) / I1;
    Sim->dynamics->getRotData(p2).angularVelocity += (r2 ^ retval.impulse) / I2;

    if (_unusedDimension != std::numeric_limits<size_t>::max()) {
      p1.getVelocity()[_unusedDimension] = 0;
      p2.getVelocity()[_unusedDimension] = 0;
      Sim->dynamics->getRotData(p1)
          .angularVelocity[(_unusedDimension + 1) % 3] = 0;
      Sim->dynamics->getRotData(p1)
          .angularVelocity[(_unusedDimension + 2) % 3] = 0;
      Sim->dynamics->getRotData(p2)
          .angularVelocity[(_unusedDimension + 1) % 3] = 0;
      Sim->dynamics->getRotData(p2)
          .angularVelocity[(_unusedDimension + 2) % 3] = 0;
    }

    return retval;
  }
  case NBHOOD_IN: {
    ICapture::add(p1, p2);
    return PairEventData(p1, p2, *Sim->species[p1], *Sim->species[p2], VIRTUAL);
  }
  case NBHOOD_OUT: {
    ICapture::remove(p1, p2);
    return PairEventData(p1, p2, *Sim->species[p1], *Sim->species[p2], VIRTUAL);
  }
  case VIRTUAL:
    return PairEventData(p1, p2, *Sim->species[p1], *Sim->species[p2], VIRTUAL);
  default:
    M_throw() << "Unknown collision type";
  }
}

void IDumbbells::outputXML(magnet::xml::XmlStream &XML) const {
  XML << magnet::xml::attr("Type") << "Dumbbells" << magnet::xml::attr("Name")
      << intName << magnet::xml::attr("Elasticity") << _e->getName();

  if (_unusedDimension != std::numeric_limits<size_t>::max())
    XML << magnet::xml::attr("UnusedDimension") << _unusedDimension;

  for (const auto &sphere : _compositeData)
    XML << magnet::xml::tag("Sphere") << magnet::xml::attr("Diameter")
        << sphere._diam->getName() << magnet::xml::tag("Offset")
        << sphere._offset / Sim->units.unitLength()
        << magnet::xml::endtag("Offset") << magnet::xml::endtag("Sphere");

  XML << range;
  ICapture::outputCaptureMap(XML);
}

size_t IDumbbells::captureTest(const Particle &p1, const Particle &p2) const {
  if (&(*(Sim->getInteraction(p1, p2))) != this)
    return false;

  const double max_dist = maxIntDist(p1, p2);
  return Sim->dynamics->sphereOverlap(p1, p2, max_dist) > 0;
}

bool IDumbbells::validateState(const Particle &p1, const Particle &p2,
                               bool textoutput) const {
  double growthrate = 0;
  if (std::dynamic_pointer_cast<DynCompression>(Sim->dynamics))
    growthrate = std::static_pointer_cast<DynCompression>(Sim->dynamics)
                     ->getGrowthRate();
  const double growthfactor = 1 + growthrate * Sim->systemTime;
  const double max_dist = maxIntDist(p1, p2) * growthfactor;

  bool has_error = false;
  double distance = Sim->BCs->getDistance(p1, p2);

  if (isCaptured(p1, p2)) {
    // Check the capture map is valid
    if (distance > max_dist) {
      if (textoutput)
        derr << "Particle " << p1.getID() << " and Particle " << p2.getID()
             << " are registered as being closer than "
             << max_dist / Sim->units.unitLength()
             << " but they're outside of this by "
             << (distance - max_dist) / Sim->units.unitLength() << std::endl;
      has_error = true;
    }

    Vector r12 = p1.getPosition() - p2.getPosition();
    Vector v12 = p1.getVelocity() - p2.getVelocity();
    Sim->BCs->applyBC(r12, v12);
    const auto &angv1 = Sim->dynamics->getRotData(p1).angularVelocity;
    const auto &angv2 = Sim->dynamics->getRotData(p2).angularVelocity;

    for (auto it1 = _compositeData.begin(); it1 != _compositeData.end(); ++it1)
      for (auto it2 = _compositeData.begin(); it2 != _compositeData.end();
           ++it2) {
        const Vector u1 =
            Sim->dynamics->getRotData(p1).orientation * it1->_offset;
        const Vector u2 =
            Sim->dynamics->getRotData(p2).orientation * it2->_offset;
        const double Diam1 = it1->_diam->getProperty(p1);
        const double Diam2 = it2->_diam->getProperty(p2);

        magnet::intersection::detail::OffcentreSpheresOverlapFunction f(
            r12, v12, angv1, angv2, u1, u2, Diam1, Diam2, max_dist,
            Sim->systemTime, growthrate, 0, 0);

        const double fval = f.eval().front();
        if (fval < 0) {
          if (textoutput)
            derr << "Comopsite particle " << p1.getID() << " sphere "
                 << it1 - _compositeData.begin() << " and Particle "
                 << p2.getID() << " sphere " << it2 - _compositeData.begin()
                 << " are overlapping by " << fval << std::endl;
          has_error = true;
        }
      }
  } else if (distance < max_dist) {
    if (textoutput)
      derr << "Particle " << p1.getID() << " and Particle " << p2.getID()
           << " are closer than " << max_dist / Sim->units.unitLength()
           << " but they've not been registered as captured, despite being at "
              "a distance of "
           << (distance - max_dist) / Sim->units.unitLength() << std::endl;
    has_error = true;
  }

  return has_error;
}
} // namespace dynamo
