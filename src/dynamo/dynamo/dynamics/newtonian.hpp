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
#include <dynamo/dynamics/dynamics.hpp>

namespace dynamo {
  /*! \brief A Dynamics which implements standard Newtonian dynamics.
    
     This Dynamics provides the dynamics of a system particles
     evolving only under interparticle Interaction(s) and Local or
     Global forces. More complex systems are available as derivations
     of this dynamics, such as a Dynamics including a constant
     gravity force (DynGravity), or a Dynamics specialized
     for multicanonical simulations (DynNewtonianMC). 
   */
  class DynNewtonian: public Dynamics
  {
  public:
    DynNewtonian(dynamo::Simulation*);

    virtual double SphereSphereInRoot(const Particle& p1, const Particle& p2, double d) const;
    virtual double SphereSphereInRoot(const Range& p1, const Range& p2, double d) const;
    virtual double SphereSphereOutRoot(const Particle& p1, const Particle& p2, double d) const;
    virtual double SphereSphereOutRoot(const Range& p1, const Range& p2, double d) const;  
    virtual double sphereOverlap(const Particle& p1, const Particle& p2, const double& d) const;
    virtual double CubeCubeInRoot(const Particle& p1, const Particle& p2, double d) const;
    virtual bool cubeOverlap(const Particle& p1, const Particle& p2, const double d) const;
    virtual void streamParticle(Particle&, const double&) const;
    virtual double getSquareCellCollision2(const Particle&, const Vector &, const Vector &) const;
    virtual int getSquareCellCollision3(const Particle&, const Vector &, const Vector &) const;
    virtual std::pair<bool,double> getPointPlateCollision(const Particle& np1, const Vector& nrw0, const Vector& nhat, const double& Delta, const double& Omega, const double& Sigma, const double& t, bool) const;
    virtual ParticleEventData runOscilatingPlate(Particle& part, const Vector& rw0, const Vector& nhat, double& delta, const double& omega0, const double& sigma, const double& mass, const double& e, double& t, bool strongPlate) const;
    virtual double getPBCSentinelTime(const Particle&, const double&) const;
    virtual PairEventData SmoothSpheresColl(const IntEvent&, const double&, const double&, const EEventType& eType) const;
    virtual bool DSMCSpheresTest(Particle&, Particle&, double&, const double&, Vector) const;
    virtual PairEventData DSMCSpheresRun(Particle&, Particle&, const double&, Vector) const;
    virtual PairEventData SphereWellEvent(const IntEvent&, const double&, const double&) const;
    virtual double getPlaneEvent(const Particle&, const Vector &, const Vector &, double) const;
    virtual std::pair<double, Dynamics::TriangleIntersectingPart> getSphereTriangleEvent(const Particle& part, const Vector & A, const Vector & B, const Vector & C, const double dist) const;
    virtual double getCylinderWallCollision(const Particle&, const Vector &, const Vector &, const double&) const;
    virtual ParticleEventData runCylinderWallCollision(Particle&, const Vector &, const Vector &, const double&) const;
    virtual ParticleEventData runPlaneEvent(Particle&, const Vector &, double, double) const;
    virtual ParticleEventData runAndersenWallCollision(Particle&, const Vector &, const double& T) const;
    virtual ParticleEventData randomGaussianEvent(Particle&, const double&, const size_t) const;
    virtual NEventData multibdyCollision(const Range&, const Range&, const double&, const EEventType&) const;
    virtual NEventData multibdyWellEvent(const Range&, const Range&, const double&, const double&, EEventType&) const;
    virtual PairEventData parallelCubeColl(const IntEvent& event, const double& e, const double& d, const EEventType& eType = CORE) const;
    virtual std::pair<bool, double> getLineLineCollision(const double length, const Particle& p1, const Particle& p2, double t_max) const;
    virtual PairEventData runLineLineCollision(const IntEvent& eevent, const double& elasticity, const double& length) const;
    virtual PairEventData RoughSpheresColl(const IntEvent& event, const double& e, const double& et, const double& d2, const EEventType& eType) const;
    virtual ParticleEventData runRoughWallCollision(Particle& part, const Vector & vNorm, const double& e, const double& et, const double& r) const;
    virtual bool getOffCenterSphereOffCenterSphereCollision(const double length, const double diameter, const Particle& p1, const Particle& p2, const double) const;
    virtual PairEventData runOffCenterSphereOffCenterSphereCollision(const IntEvent& eevent, const double& elasticity, const double& length, const double& diameter) const;

  protected:
    virtual void outputXML(magnet::xml::XmlStream&) const;

    mutable long double lastAbsoluteClock;
    mutable unsigned int lastCollParticle1;
    mutable unsigned int lastCollParticle2;
  };
}
