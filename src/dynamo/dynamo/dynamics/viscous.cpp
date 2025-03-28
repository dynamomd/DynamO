/*  dynamo:- Event driven molecular dynamics simulator
    http://www.dynamomd.org
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>
    Copyright (C) 2011  Sebastian Gonzalez

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
#include <dynamo/dynamics/viscous.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/species/species.hpp>
#include <dynamo/units/units.hpp>
#include <magnet/intersection/line_line.hpp>
#include <magnet/intersection/overlapfuncs/oscillatingplate.hpp>
#include <magnet/intersection/ray_cube.hpp>
#include <magnet/intersection/ray_plane.hpp>
#include <magnet/intersection/ray_rod.hpp>
#include <magnet/intersection/ray_sphere.hpp>
#include <magnet/intersection/ray_triangle.hpp>
#include <magnet/math/matrix.hpp>
#include <magnet/overlap/point_cube.hpp>
#include <magnet/overlap/point_prism.hpp>
#include <magnet/xmlwriter.hpp>

#include <magnet/intersection/polynomial.hpp>

namespace dynamo {
DynViscous::DynViscous(dynamo::Simulation *tmp, const magnet::xml::Node &XML)
    : DynNewtonian(tmp), _g({0, -1, 0}) {
  _g << XML.getNode("g");
  _g *= Sim->units.unitAcceleration();
  _gamma =
      XML.getAttribute("gamma").as<double>() * Sim->units.unitAcceleration();
}

void DynViscous::outputXML(magnet::xml::XmlStream &XML) const {
  XML << magnet::xml::attr("Type") << "Viscous";
  XML << magnet::xml::attr("gamma") << _gamma * Sim->units.unitTime();
  XML << magnet::xml::tag("g") << _g / Sim->units.unitAcceleration()
      << magnet::xml::endtag("g");
}

void DynViscous::streamParticle(Particle &particle, const double &dt) const {
  const double m = Sim->species[particle]->getMass(particle.getID());

  particle.getPosition() +=
      -_g * dt / (_gamma * m) + (1 - std::exp(-_gamma * dt)) *
                                    (particle.getVelocity() + _g / _gamma) /
                                    (_gamma * m);

  particle.getVelocity() =
      -_g / _gamma +
      (particle.getVelocity() + _g / _gamma) * std::exp(-_gamma * dt);

  // This part is also incorrect, but left here for future work
  if (hasOrientationData()) {
    orientationData[particle.getID()].orientation =
        Quaternion::fromRotationAxis(
            orientationData[particle.getID()].angularVelocity * dt) *
        orientationData[particle.getID()].orientation;
    orientationData[particle.getID()].orientation.normalise();
  }
}

double DynViscous::SphereSphereInRoot(const Particle &p1, const Particle &p2,
                                      double sigma) const {
  Vector r12 = p1.getPosition() - p2.getPosition();
  Sim->BCs->applyBC(r12);
  const Vector X = r12;
  const double m1 = Sim->species[p1]->getMass(p1.getID());
  const double m2 = Sim->species[p2]->getMass(p2.getID());
  const Vector V = (p1.getVelocity() / m1) - (p2.getVelocity() / m2);

  if (m1 != m2)
    M_throw()
        << "Not implemented asymmetric particle masses for viscous dynamics";

  const double c = X.nrm2() - sigma * sigma;
  const double b = 2 * (V | X) / _gamma;
  const double a = V.nrm2() / (_gamma * _gamma);

  // As the transform from y=e^{-\gamma t} to t is a monotonic
  // transform we can perform the stable algorithm in y, but we need
  // to ensure that t=0 corresponds to y=0 and t>0 corresponds to y>0
  //  Thus the appropriate transformation is y=1-e^{-\gamma\,t}
  magnet::intersection::detail::PolynomialFunction<2> f(c, b, 2 * a);
  const double y = magnet::intersection::detail::nextEvent(f);

  // If y>1 then there's no real roots
  if (y > 1)
    return HUGE_VAL;

  const double t = -std::log(1 - y) / _gamma;
  return t;
}

PairEventData DynViscous::SmoothSpheresColl(Event &event, const double &e,
                                            const double &,
                                            const EEventType &eType) const {
  Particle &particle1 = Sim->particles[event._particle1ID];
  Particle &particle2 = Sim->particles[event._particle2ID];
  updateParticlePair(particle1, particle2);
  PairEventData retVal(particle1, particle2, *Sim->species[particle1],
                       *Sim->species[particle2], eType);

  Sim->BCs->applyBC(retVal.rij, retVal.vijold);

  double p1Mass = Sim->species[retVal.particle1_.getSpeciesID()]->getMass(
      particle1.getID());
  double p2Mass = Sim->species[retVal.particle2_.getSpeciesID()]->getMass(
      particle2.getID());

  retVal.rvdot = retVal.rij | retVal.vijold;

  double mu = 1.0 / ((1.0 / p1Mass) + (1.0 / p2Mass));

  // If both particles have infinite mass, we need to modify the
  // masses (and mu) to allow collisions.
  bool infinite_masses = (p1Mass == std::numeric_limits<float>::infinity()) &&
                         (p2Mass == std::numeric_limits<float>::infinity());
  if (infinite_masses) {
    p1Mass = p2Mass = 1;
    mu = 0.5;
  }

  retVal.impulse =
      retVal.rij * ((1.0 + e) * mu * retVal.rvdot / retVal.rij.nrm2());

  bool neededFix = false;
  while (true) {
    const Vector v1n = particle1.getVelocity() - retVal.impulse / p1Mass;
    const Vector v2n = particle2.getVelocity() + retVal.impulse / p2Mass;
    Vector rijn = particle1.getPosition() - particle2.getPosition();
    Vector vijn = v1n - v2n;
    Sim->BCs->applyBC(rijn, vijn);
    if (((retVal.rvdot < 0) && ((rijn | vijn) > 0)) ||
        ((retVal.rvdot > 0) && ((rijn | vijn) < 0))) {
      particle1.getVelocity() = v1n;
      particle2.getVelocity() = v2n;
      break;
    }
    neededFix = true;
    retVal.impulse *= 2;
  }
  if (neededFix)
    retVal.impulse *= 8;

  retVal.impulse *= !infinite_masses;

  lastCollParticle1 = particle1.getID();
  lastCollParticle2 = particle2.getID();
  lastAbsoluteClock = Sim->systemTime;
  return retVal;
}

double DynViscous::getPBCSentinelTime(const Particle &part,
                                      const double &lMax) const {
  // This is bad for low densities!
  return HUGE_VAL;
}
} // namespace dynamo
