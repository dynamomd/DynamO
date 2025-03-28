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

#pragma once
#include <dynamo/2particleEventData.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/dynamics/newtonian.hpp>

namespace dynamo {
/*! \brief A Dynamics which implements "viscous" Newtonian dynamics.

   This Dynamics provides the dynamics of a system of particles
   damped by a viscous term.
 */
class DynViscous : public DynNewtonian {
public:
  DynViscous(dynamo::Simulation *, const magnet::xml::Node &);
  virtual double SphereSphereInRoot(const Particle &p1, const Particle &p2,
                                    double d) const;
  virtual void streamParticle(Particle &, const double &) const;
  virtual double getPBCSentinelTime(const Particle &, const double &) const;
  virtual PairEventData SmoothSpheresColl(Event &, const double &,
                                          const double &,
                                          const EEventType &eType) const;

  Vector _g;
  double _gamma;

  // Inherited
  //  virtual bool cubeOverlap(const Particle& p1, const Particle& p2, const
  //  double d) const; virtual double sphereOverlap(const Particle& p1, const
  //  Particle& p2, const double& d) const; virtual PairEventData
  //  RoughSpheresColl(Event& event, const double& e, const double& et, const
  //  double& d1, const double& d2, const EEventType& eType) const { M_throw()
  //  << "Not implemented";} virtual PairEventData SmoothSpheresColl(Event&,
  //  const double&, const double&, const EEventType& eType) const;

  // Not implemented
  virtual double SphereSphereInRoot(const IDRange &p1, const IDRange &p2,
                                    double d) const {
    M_throw() << "Not implemented";
  }
  virtual double SphereSphereOutRoot(const Particle &p1, const Particle &p2,
                                     double d) const {
    M_throw() << "Not implemented";
  }
  virtual double SphereSphereOutRoot(const IDRange &p1, const IDRange &p2,
                                     double d) const {
    M_throw() << "Not implemented";
  }
  virtual double CubeCubeInRoot(const Particle &p1, const Particle &p2,
                                double d) const {
    M_throw() << "Not implemented";
  }
  virtual double getSquareCellCollision2(const Particle &, const Vector &,
                                         const Vector &) const {
    M_throw() << "Not implemented";
  }
  virtual int getSquareCellCollision3(const Particle &, const Vector &,
                                      const Vector &) const {
    M_throw() << "Not implemented";
  }
  virtual std::pair<bool, double>
  getPointPlateCollision(const Particle &np1, const Vector &nrw0,
                         const Vector &nhat, const double &Delta,
                         const double &Omega, const double &Sigma,
                         const double &t, bool) const {
    M_throw() << "Not implemented";
  }
  virtual ParticleEventData
  runOscilatingPlate(Particle &part, const Vector &rw0, const Vector &nhat,
                     double &delta, const double &omega0, const double &sigma,
                     const double &mass, const double &e, double &t,
                     bool strongPlate) const {
    M_throw() << "Not implemented";
  }
  virtual bool DSMCSpheresTest(Particle &, Particle &, double &, const double &,
                               Vector) const {
    M_throw() << "Not implemented";
  }
  virtual PairEventData DSMCSpheresRun(Particle &, Particle &, const double &,
                                       Vector) const {
    M_throw() << "Not implemented";
  }
  virtual PairEventData SphereWellEvent(Event &, const double &, const double &,
                                        size_t) const {
    M_throw() << "Not implemented";
  }
  virtual double getPlaneEvent(const Particle &, const Vector &, const Vector &,
                               double) const {
    M_throw() << "Not implemented";
  }
  virtual ParticleEventData runPlaneEvent(Particle &, const Vector &,
                                          const double, const double) const {
    M_throw() << "Not implemented";
  }
  virtual std::pair<double, Dynamics::TriangleIntersectingPart>
  getSphereTriangleEvent(const Particle &part, const Vector &A, const Vector &B,
                         const Vector &C, const double dist) const {
    M_throw() << "Not implemented";
  }
  virtual double getCylinderWallCollision(const Particle &, const Vector &,
                                          const Vector &,
                                          const double &) const {
    M_throw() << "Not implemented";
  }
  virtual ParticleEventData runCylinderWallCollision(Particle &, const Vector &,
                                                     const Vector &,
                                                     const double &) const {
    M_throw() << "Not implemented";
  }
  virtual ParticleEventData runAndersenWallCollision(Particle &, const Vector &,
                                                     const double &T,
                                                     const double d,
                                                     const double slip) const {
    M_throw() << "Not implemented";
  }
  virtual ParticleEventData randomGaussianEvent(Particle &, const double &,
                                                const size_t) const {
    M_throw() << "Not implemented";
  }
  virtual NEventData multibdyCollision(const IDRange &, const IDRange &,
                                       const double &,
                                       const EEventType &) const {
    M_throw() << "Not implemented";
  }
  virtual NEventData multibdyWellEvent(const IDRange &, const IDRange &,
                                       const double &, const double &,
                                       EEventType &) const {
    M_throw() << "Not implemented";
  }
  virtual PairEventData parallelCubeColl(Event &event, const double &e,
                                         const double &d,
                                         const EEventType &eType = CORE) const {
    M_throw() << "Not implemented";
  }
  virtual std::pair<bool, double> getLineLineCollision(const double length,
                                                       const Particle &p1,
                                                       const Particle &p2,
                                                       double t_max) const {
    M_throw() << "Not implemented";
  }
  virtual PairEventData runLineLineCollision(Event &eevent,
                                             const double &elasticity,
                                             const double &length) const {
    M_throw() << "Not implemented";
  }
  virtual ParticleEventData
  runRoughWallCollision(Particle &part, const Vector &vNorm, const double &e,
                        const double &et, const double &r) const {
    M_throw() << "Not implemented";
  }

protected:
  virtual void outputXML(magnet::xml::XmlStream &) const;

  mutable long double lastAbsoluteClock;
  mutable unsigned int lastCollParticle1;
  mutable unsigned int lastCollParticle2;
};
} // namespace dynamo
