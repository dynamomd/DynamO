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

#include <dynamo/dynamics/viscous.hpp>
#include <dynamo/2particleEventData.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/BC/BC.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/species/species.hpp>
#include <dynamo/units/units.hpp>
#include <magnet/overlap/point_prism.hpp>
#include <magnet/overlap/point_cube.hpp>
#include <magnet/intersection/ray_triangle.hpp>
#include <magnet/intersection/ray_rod.hpp>
#include <magnet/intersection/ray_sphere.hpp>
#include <magnet/intersection/ray_plane.hpp>
#include <magnet/intersection/ray_cube.hpp>
#include <magnet/intersection/line_line.hpp>
#include <magnet/intersection/overlapfuncs/oscillatingplate.hpp>
#include <magnet/math/matrix.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
  bool 
  DynViscous::cubeOverlap(const Particle& p1, const Particle& p2, 
			    const double d) const
  {
    Vector r12 = p1.getPosition() - p2.getPosition();
    Sim->BCs->applyBC(r12);
    return magnet::overlap::point_cube(r12, 2 * Vector{d, d, d});
  }

  double
  DynViscous::SphereSphereInRoot(const Particle& p1, const Particle& p2, double d) const
  {
    Vector r12 = p1.getPosition() - p2.getPosition();
    Vector v12 = p1.getVelocity() - p2.getVelocity();
    Sim->BCs->applyBC(r12, v12);
    return magnet::intersection::ray_sphere(r12, v12, d);
  }
  
  DynViscous::DynViscous(dynamo::Simulation* tmp, const magnet::xml::Node& xml):
    Dynamics(tmp), lastAbsoluteClock(-1), lastCollParticle1(0), lastCollParticle2(0)  
  {
  }

  void
  DynViscous::streamParticle(Particle &particle, const double &dt) const
  {
    particle.getPosition() += particle.getVelocity() * dt;

    if (hasOrientationData())
      {
	orientationData[particle.getID()].orientation = Quaternion::fromRotationAxis(orientationData[particle.getID()].angularVelocity * dt)
	  * orientationData[particle.getID()].orientation ;
	orientationData[particle.getID()].orientation.normalise();
      }
  }

  double 
  DynViscous::getPlaneEvent(const Particle& part, const Vector& wallLoc, const Vector& wallNorm, double diameter) const
  {
    Vector rij = part.getPosition() - wallLoc, vel = part.getVelocity();
    Sim->BCs->applyBC(rij, vel);

    
    return magnet::intersection::ray_plane(rij, vel, wallNorm, diameter);
  }

  ParticleEventData 
  DynViscous::runPlaneEvent(Particle& part, const Vector& vNorm, const double e, const double diameter) const
  {
    updateParticle(part);

    ParticleEventData retVal(part, *Sim->species(part), WALL);
  
    part.getVelocity() -= (1+e) * (vNorm | part.getVelocity()) * vNorm;
  
    return retVal; 
  }

  double
  DynViscous::getSquareCellCollision2(const Particle& part, const Vector & origin, const Vector & width) const
  {
    Vector  rpos(part.getPosition() - origin);
    Vector  vel(part.getVelocity());
    Sim->BCs->applyBC(rpos, vel);
  
    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      if ((vel[iDim] == 0) && (std::signbit(vel[iDim])))
	vel[iDim] = 0;

    double retVal;
    if (vel[0] < 0)
      retVal = -rpos[0] / vel[0];
    else
      retVal = (width[0]-rpos[0]) / vel[0];

    for (size_t iDim = 1; iDim < NDIM; ++iDim)
      {
	double tmpdt((vel[iDim] < 0)
		     ? -rpos[iDim]/vel[iDim] 
		     : (width[iDim]-rpos[iDim]) / vel[iDim]);
      
	if (tmpdt < retVal)
	  retVal = tmpdt;
      }
  
    return retVal;
  }

  int
  DynViscous::getSquareCellCollision3(const Particle& part, const Vector & origin, const Vector & width) const
  {
    Vector rpos(part.getPosition() - origin);
    Vector vel(part.getVelocity());
    Sim->BCs->applyBC(rpos, vel);

    int retVal(0);
    double time(std::numeric_limits<float>::infinity());
  
    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      if ((vel[iDim] == 0) && (std::signbit(vel[iDim])))
	vel[iDim] = 0;

    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      {
	double tmpdt = ((vel[iDim] < 0) ? -rpos[iDim]/vel[iDim] : (width[iDim]-rpos[iDim]) / vel[iDim]);

	if (tmpdt < time)
	  {
	    time = tmpdt;
	    retVal = (vel[iDim] < 0) ? -(iDim+1) : (iDim+1);
	  }
      }

    if (((retVal < 0) && (vel[abs(retVal)-1] > 0)) || ((retVal > 0) && (vel[abs(retVal)-1] < 0)))
      M_throw() << "Found an error! retVal " << retVal
		<< " vel is " << vel[abs(retVal)-1];

    return retVal;
  }

  PairEventData 
  DynViscous::SmoothSpheresColl(Event& event, const double& e, const double&, const EEventType& eType) const
  {
    Particle& particle1 = Sim->particles[event._particle1ID];
    Particle& particle2 = Sim->particles[event._particle2ID];
    updateParticlePair(particle1, particle2);  
    PairEventData retVal(particle1, particle2, *Sim->species(particle1), *Sim->species(particle2), eType);

    Sim->BCs->applyBC(retVal.rij, retVal.vijold);

    double p1Mass = Sim->species[retVal.particle1_.getSpeciesID()]->getMass(particle1.getID()); 
    double p2Mass = Sim->species[retVal.particle2_.getSpeciesID()]->getMass(particle2.getID());
 
    retVal.rvdot = retVal.rij | retVal.vijold;

    double mu = 1.0 / ((1.0 / p1Mass) + (1.0 / p2Mass));

    //If both particles have infinite mass, we need to modify the
    //masses (and mu) to allow collisions.
    bool infinite_masses = (p1Mass == std::numeric_limits<float>::infinity()) && (p2Mass == std::numeric_limits<float>::infinity());
    if (infinite_masses)
      {
	p1Mass = p2Mass = 1;
	mu = 0.5;
      }

    retVal.impulse = retVal.rij * ((1.0 + e) * mu * retVal.rvdot / retVal.rij.nrm2());  
    particle1.getVelocity() -= retVal.impulse / p1Mass;
    particle2.getVelocity() += retVal.impulse / p2Mass;
    retVal.impulse *= !infinite_masses;

    lastCollParticle1 = particle1.getID();
    lastCollParticle2 = particle2.getID();
    lastAbsoluteClock = Sim->systemTime;
    return retVal;
  }

  void 
  DynViscous::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "Viscous";
  }

  double 
  DynViscous::getPBCSentinelTime(const Particle& part, const double& lMax) const
  {
#ifdef DYNAMO_DEBUG
    if (!isUpToDate(part))
      M_throw() << "Particle is not up to date";
#endif

    Vector pos(part.getPosition()), vel(part.getVelocity());

    Sim->BCs->applyBC(pos, vel);

    double retval = std::numeric_limits<float>::infinity();

    for (size_t i(0); i < NDIM; ++i)
      if (vel[i] != 0)
	{
	  double tmp = (0.5 * (0.5 * Sim->primaryCellSize[i] - lMax)) / std::abs(vel[i]);
	  
	  if (tmp < retval)
	    retval = tmp;
	}

    return retval;
  }

  double 
  DynViscous::sphereOverlap(const Particle& p1, const Particle& p2, const double& d) const
  {
    Vector r12 = p1.getPosition() - p2.getPosition();
    Sim->BCs->applyBC(r12);

    return std::max(d  - std::sqrt(r12 | r12), 0.0);
  }

  PairEventData 
  DynViscous::RoughSpheresColl(Event& event, const double& e, const double& et, const double& d1, const double& d2, const EEventType& eType) const
  {
    if (!hasOrientationData())
      M_throw() << "Cannot use tangential coefficients of inelasticity without orientational data/species";

    Particle& particle1 = Sim->particles[event._particle1ID];
    Particle& particle2 = Sim->particles[event._particle2ID];

    updateParticlePair(particle1, particle2);  
    PairEventData retVal(particle1, particle2, *Sim->species(particle1), *Sim->species(particle2), eType);
    
    Sim->BCs->applyBC(retVal.rij, retVal.vijold);
  
    double p1Mass = Sim->species[retVal.particle1_.getSpeciesID()]->getMass(particle1.getID());
    double p2Mass = Sim->species[retVal.particle2_.getSpeciesID()]->getMass(particle2.getID());

    retVal.rvdot = (retVal.rij | retVal.vijold);

    const Vector rijhat = retVal.rij / retVal.rij.nrm();
    const Vector gij = retVal.vijold - ((0.5 * d1 * orientationData[particle1.getID()].angularVelocity + 0.5 * d2 * orientationData[particle2.getID()].angularVelocity) ^ rijhat);
    const Vector rcrossgij = rijhat ^ gij;
    const double rdotgij = rijhat | gij;

    double mu = 1.0 / ((1.0 / p1Mass) + (1.0 / p2Mass));

    double I = 2.0/5.0;
    
    bool infinite_masses = (p1Mass == std::numeric_limits<float>::infinity()) && (p2Mass == std::numeric_limits<float>::infinity());
    if (infinite_masses)
      {
	p1Mass = p2Mass = 1;
	mu = 0.5;
      }
    
    retVal.impulse = mu * ((1+e) * rijhat * rdotgij + ((et - 1) / (1 + 1.0 / I)) * (rijhat ^ (rcrossgij)));
    particle1.getVelocity() -= retVal.impulse / p1Mass;
    particle2.getVelocity() += retVal.impulse / p2Mass;

    retVal.impulse *= !infinite_masses;

    Vector angularVchange = (mu * (1-et) / (1 + I)) * rcrossgij;
 
    orientationData[particle1.getID()].angularVelocity += angularVchange / (p1Mass * d1 * 0.5);
    orientationData[particle2.getID()].angularVelocity += angularVchange / (p2Mass * d2 * 0.5);
    return retVal;
  }
}
