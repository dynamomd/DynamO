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

#include <dynamo/dynamics/compression.hpp>
#include <dynamo/interactions/intEvent.hpp>
#include <dynamo/2particleEventData.hpp>

#include <dynamo/BC/BC.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/species/species.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <magnet/intersection/ray_plane.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
  DynCompression::DynCompression(dynamo::Simulation* tmp, double GR):
    DynNewtonian(tmp),
    growthRate(GR) {}

  double 
  DynCompression::SphereSphereInRoot(const Particle& p1, const Particle& p2, double d) const
  {
    Vector r12 = p1.getPosition() - p2.getPosition();
    Vector v12 = p1.getVelocity() - p2.getVelocity();
    Sim->BCs->applyBC(r12, v12);
    double rvdot = r12 | v12;
    double r2 = r12 | r12;
    double v2 = v12 | v12;
    double d2 = d * d;
    
    double b = rvdot - d2 
      * (growthRate * growthRate * Sim->systemTime + growthRate);
  
    if (b >= 0.0) return HUGE_VAL;


    double a = v2 - growthRate * growthRate * d2;
    double c = r2 - d2 * (1.0 + growthRate * Sim->systemTime 
			  * (2.0 + growthRate * Sim->systemTime));
    double arg = (b * b) - a * c;
    
    if (arg < 0.0) return HUGE_VAL;

    return  c / (sqrt(arg) - b);
  }
  
  double 
  DynCompression::getPlaneEvent(const Particle& part, const Vector & origin, const Vector & norm, double diameter) const
  {
    Vector rij = part.getPosition() - origin,
      vij = part.getVelocity() - norm * diameter * growthRate;

    Sim->BCs->applyBC(rij, vij);

    return magnet::intersection::ray_plane(rij, vij, norm, diameter * (1.0 + growthRate * Sim->systemTime));
  }

  ParticleEventData 
  DynCompression::runPlaneEvent(Particle& part, const Vector& vNorm, const double e, const double diameter) const
  {
    updateParticle(part);

    ParticleEventData retVal(part, *Sim->species[part], WALL);
  
    Vector vij = part.getVelocity() - vNorm * diameter * growthRate;

    part.getVelocity() -= (1+e) * (vNorm | vij) * vNorm;
  
    retVal.setDeltaKE(0.5 * Sim->species[retVal.getSpeciesID()]->getMass(part.getID())
		      * (part.getVelocity().nrm2() - retVal.getOldVel().nrm2()));
    return retVal; 
  }

  double 
  DynCompression::SphereSphereOutRoot(const Particle& p1, const Particle& p2, double d) const
  {
    Vector r12 = p1.getPosition() - p2.getPosition();
    Vector v12 = p1.getVelocity() - p2.getVelocity();
    Sim->BCs->applyBC(r12, v12);
    double rvdot = r12 | v12;
    double r2 = r12 | r12;
    double v2 = v12 | v12;
    double d2 = d * d;

    double a = v2 - growthRate * growthRate * d2;
    double b = rvdot - d2 * (growthRate * growthRate 
				 * Sim->systemTime + growthRate);
    double c = d2 * (1.0 + growthRate * Sim->systemTime 
		     * (2.0 + growthRate * Sim->systemTime)) - r2;

    double arg = (b * b) + a * c;
  
    if (arg > 0.0)
      if (a > 0.0)
	{
	  if (b < 0)
	    return (sqrt(arg) - b) / a;
	  
	  return c / (sqrt(arg) + b);
	}
    
    return HUGE_VAL;
  }

  double
  DynCompression::sphereOverlap(const Particle& p1, const Particle& p2,
			      const double& d) const
  {
    Vector r12 = p1.getPosition() - p2.getPosition();
    Sim->BCs->applyBC(r12);

    double currd2 = d * d * (1 + 2.0 * Sim->systemTime * growthRate 
			     + pow(Sim->systemTime * growthRate, 2));
    return std::sqrt(std::max(currd2 - (r12 | r12), 0.0));
  }


  PairEventData 
  DynCompression::SmoothSpheresColl(const IntEvent& event, const double& e, const double& d2, const EEventType& eType) const
  {
    Particle& particle1 = Sim->particles[event.getParticle1ID()];
    Particle& particle2 = Sim->particles[event.getParticle2ID()];

    updateParticlePair(particle1, particle2);  

    PairEventData retVal(particle1, particle2,
			 *Sim->species[particle1],
			 *Sim->species[particle2],
			 eType);

    Sim->BCs->applyBC(retVal.rij, retVal.vijold);
    
    double p1Mass = Sim->species[retVal.particle1_.getSpeciesID()]->getMass(particle1.getID()); 
    double p2Mass = Sim->species[retVal.particle2_.getSpeciesID()]->getMass(particle2.getID()); 
    double r2 = retVal.rij.nrm2();
  
    retVal.rvdot = (retVal.rij | retVal.vijold);

    //Treat special cases if one particle has infinite mass
    if ((p1Mass == 0) && (p2Mass != 0))
      {
	retVal.impulse = p2Mass * retVal.rij * ((1.0 + e) * (retVal.rvdot - growthRate * sqrt(d2 * r2)) / retVal.rij.nrm2());  
	particle2.getVelocity() += retVal.impulse / p2Mass;
      }
    else if ((p2Mass == 0) && (p1Mass != 0))
      {
	retVal.impulse = p1Mass * retVal.rij * ((1.0 + e) * (retVal.rvdot - growthRate * sqrt(d2 * r2)) / retVal.rij.nrm2());  
	particle1.getVelocity() -= retVal.impulse / p1Mass;
      }
    else
      {
	bool isInfInf = ((p1Mass == 0.0) && (p2Mass == 0.0));

	//If both particles have infinite mass we just collide them as identical masses
	double mu = isInfInf ? 0.5 : p1Mass * p2Mass / (p1Mass + p2Mass);

	retVal.impulse = retVal.rij * ((1.0 + e) * mu * (retVal.rvdot - growthRate * sqrt(d2 * r2)) / retVal.rij.nrm2());  

	particle1.getVelocity() -= retVal.impulse / (p1Mass + isInfInf);
	particle2.getVelocity() += retVal.impulse / (p2Mass + isInfInf);

	//If both particles have infinite mass we pretend no momentum was transferred
	retVal.impulse *= !isInfInf;
      }


    retVal.particle1_.setDeltaKE(0.5 * p1Mass
				 * (particle1.getVelocity().nrm2() 
				    - retVal.particle1_.getOldVel().nrm2()));
  
    retVal.particle2_.setDeltaKE(0.5 * p2Mass
				 * (particle2.getVelocity().nrm2()
				    - retVal.particle2_.getOldVel().nrm2()));
  
    return retVal;
  }

  PairEventData 
  DynCompression::SphereWellEvent(const IntEvent& event, const double& deltaKE, const double& d2, size_t) const
  {
    Particle& particle1 = Sim->particles[event.getParticle1ID()];
    Particle& particle2 = Sim->particles[event.getParticle2ID()];

    updateParticlePair(particle1, particle2);  
  
    PairEventData retVal(particle1, particle2,
			 *Sim->species[particle1],
			 *Sim->species[particle2],
			 event.getType());
    
    Sim->BCs->applyBC(retVal.rij, retVal.vijold);
    
    double p1Mass = Sim->species[retVal.particle1_.getSpeciesID()]->getMass(particle1.getID());
    double p2Mass = Sim->species[retVal.particle2_.getSpeciesID()]->getMass(particle2.getID());
    double mu = p1Mass * p2Mass / (p1Mass + p2Mass);  
    Vector  urij = retVal.rij / retVal.rij.nrm();

    retVal.rvdot = (urij | retVal.vijold);

    double sqrtArg = std::pow(retVal.rvdot - growthRate * sqrt(d2), 2)  + (2.0 * deltaKE / mu);
    if ((deltaKE < 0) && (sqrtArg < 0))
      {
	event.setType(BOUNCE);
	retVal.setType(BOUNCE);

	retVal.impulse = urij * (2.0 * mu * (retVal.rvdot - growthRate * sqrt(d2)));
      }
    else if (deltaKE==0)
      retVal.impulse = Vector(0,0,0);
    else
      {	  
	retVal.particle1_.setDeltaU(-0.5 * deltaKE);
	retVal.particle2_.setDeltaU(-0.5 * deltaKE);	  
      
	if (retVal.rvdot < 0)
	  retVal.impulse = urij 
	    * (2.0 * deltaKE / (growthRate * sqrt(d2) + std::sqrt(sqrtArg) - retVal.rvdot ));
	else
	  retVal.impulse = urij 
	    * (2.0 * deltaKE / (growthRate * sqrt(d2) - std::sqrt(sqrtArg) - retVal.rvdot ));
	;
      }

    retVal.rvdot *= retVal.rij.nrm();
  
#ifdef DYNAMO_DEBUG
    if (boost::math::isnan(retVal.impulse[0]))
      M_throw() << "A nan dp has ocurred"
		<< "\ndeltaKE = " << deltaKE
		<< "\ngrowthRate = " << growthRate
		<< "\nd2 = " << d2
		<< "\nsqrtArg = " << sqrtArg
		<< "\nrvdot = " << retVal.rvdot
		<< "\nArg " << (growthRate * sqrt(d2) - std::sqrt(sqrtArg) - retVal.rvdot)
	;
#endif
  
    //This function must edit particles so it overrides the const!
    particle1.getVelocity() -= retVal.impulse / p1Mass;
    particle2.getVelocity() += retVal.impulse / p2Mass;
  
    retVal.particle1_.setDeltaKE(0.5 * p1Mass
				 * (particle1.getVelocity().nrm2() 
				    - retVal.particle1_.getOldVel().nrm2()));
  
    retVal.particle2_.setDeltaKE(0.5 * p2Mass
				 * (particle2.getVelocity().nrm2() 
				    - retVal.particle2_.getOldVel().nrm2()));
  
    return retVal;
  }

  void 
  DynCompression::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") 
	<< "Compression";
  }

  double 
  DynCompression::getPBCSentinelTime(const Particle& part,
				   const double& lMax) const
  {
#ifdef DYNAMO_DEBUG
    if (!isUpToDate(part))
      M_throw() << "Particle is not up to date";
#endif

    Vector  pos(part.getPosition()), vel(part.getVelocity());

    Sim->BCs->applyBC(pos, vel);

    double retval = (0.5 * Sim->primaryCellSize[0] - lMax) / (fabs(vel[0]) + lMax * growthRate);

    for (size_t i(1); i < NDIM; ++i)
      {
	double tmp = (0.5 * Sim->primaryCellSize[i] - lMax) / (fabs(vel[0]) + lMax * growthRate);
      
	if (tmp < retval)
	  retval = tmp;
      }

    return retval;
  }

  PairEventData 
  DynCompression::parallelCubeColl(const IntEvent&,
				 const double&, const double&,
				 const EEventType&) const
  { M_throw() << "Not Implemented"; }


  ParticleEventData 
  DynCompression::runAndersenWallCollision(Particle& part, const Vector & vNorm, const double& sqrtT, const double d) const
  {  
    updateParticle(part);

    if (hasOrientationData())
      M_throw() << "Need to implement thermostating of the rotational degrees"
	" of freedom";

    //This gives a completely new random unit vector with a properly
    //distributed Normal component. See Granular Simulation Book
    ParticleEventData tmpDat(part, *Sim->species[part], WALL);
 
    double mass = Sim->species[tmpDat.getSpeciesID()]->getMass(part.getID());

    for (size_t iDim = 0; iDim < NDIM; iDim++)
      part.getVelocity()[iDim] = Sim->normal_sampler() * sqrtT / std::sqrt(mass);
  
    part.getVelocity() 
      //This first line adds a component in the direction of the normal
      += vNorm * (sqrtT * sqrt(-2.0*log(1.0-Sim->uniform_sampler()) / mass)
		  //This removes the original normal component
		  -(part.getVelocity() | vNorm)
		  //This adds on the velocity of the wall
		  + d * growthRate)
      ;

    tmpDat.setDeltaKE(0.5 * mass * (part.getVelocity().nrm2() - tmpDat.getOldVel().nrm2()));
  
    return tmpDat; 
  }
}
