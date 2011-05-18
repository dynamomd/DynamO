/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
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

#include "NewtonL.hpp"
#include "../interactions/intEvent.hpp"
#include "../2particleEventData.hpp"
#include "../NparticleEventData.hpp"
#include "../dynamics.hpp"
#include "../BC/BC.hpp"
#include "../../base/is_simdata.hpp"
#include "../species/species.hpp"
#include "../../schedulers/sorters/datastruct.hpp"
#include "shapes/frenkelroot.hpp"
#include "shapes/oscillatingplate.hpp"
#include "shapes/lines.hpp"
#include "shapes/dumbbells.hpp"
#include "../units/units.hpp"
#include <boost/math/special_functions/fpclassify.hpp>
#include <magnet/math/matrix.hpp>
#include <magnet/xmlwriter.hpp>

bool 
LNewtonian::CubeCubeInRoot(CPDData& dat, const double& d) const
{
  //To be approaching, the largest dimension of rij must be shrinking
  
  size_t largedim(0);
  for (size_t iDim(1); iDim < NDIM; ++iDim)
    if (fabs(dat.rij[iDim]) > fabs(dat.rij[largedim])) largedim = iDim;
    
  if (dat.rij[largedim] * dat.vij[largedim] >= 0) return false;

  double tInMax(-HUGE_VAL), tOutMin(HUGE_VAL);
  
  for (size_t iDim(0); iDim < NDIM; ++iDim)
    {
      double tmptime1 = -(dat.rij[iDim] + d) / dat.vij[iDim];
      double tmptime2 = -(dat.rij[iDim] - d) / dat.vij[iDim];
      
      if (tmptime1 < tmptime2)
	{
	  if (tmptime1 > tInMax) tInMax = tmptime1;
	  if (tmptime2 < tOutMin) tOutMin = tmptime2;
	}
      else
	{
	  if (tmptime2 > tInMax) tInMax = tmptime2;
	  if (tmptime1 < tOutMin) tOutMin = tmptime1;
	}
    }
  
  if (tInMax >= tOutMin) return false;

  dat.dt = tInMax;
  return true;
}

bool 
LNewtonian::cubeOverlap(const CPDData& dat, const double& d) const
{
  for (size_t iDim(0); iDim < NDIM; ++iDim)
    if (fabs(dat.rij[iDim]) > d) return false;
  
  return true;
}

bool 
LNewtonian::SphereSphereInRoot(CPDData& dat, const double& d2, bool, bool) const
{
  if (dat.rvdot < 0)
    {
      double arg = dat.rvdot * dat.rvdot - dat.v2 * (dat.r2 - d2);

      if (arg > 0)
	{
	  //This is the more numerically stable form of the quadratic
	  //formula
	  dat.dt = (d2 - dat.r2) / (dat.rvdot-sqrt(arg));

#ifdef DYNAMO_DEBUG
	  if (boost::math::isnan(dat.dt))
	    M_throw() << "dat.dt is nan";
#endif
	  return true;
	}
      else 
	return false;
    }
  else
    return false;
}
  
bool 
LNewtonian::SphereSphereOutRoot(CPDData& dat, const double& d2, bool, bool) const
{
  dat.dt = (sqrt(dat.rvdot * dat.rvdot 
		 - dat.v2 * (dat.r2 - d2)) - dat.rvdot) / dat.v2;

  if (boost::math::isnan(dat.dt))
    {//The nan occurs if the spheres aren't moving apart
      dat.dt = HUGE_VAL;
      return false;
    }
  else
    return true;
}

bool 
LNewtonian::sphereOverlap(const CPDData& dat, const double& d2) const
{
  return (dat.r2 - d2) < 0.0;
}

ParticleEventData 
LNewtonian::randomGaussianEvent(const Particle& part, const double& sqrtT) const
{
  //See http://mathworld.wolfram.com/SpherePointPicking.html

  if (hasOrientationData())
    M_throw() << "Need to implement thermostating of the rotational degrees"
      " of freedom";  

  //Ensure the particle is free streamed first
  updateParticle(part);

  //Collect the precoll data
  ParticleEventData tmpDat(part, Sim->dynamics.getSpecies(part), GAUSSIAN);

  double mass = tmpDat.getSpecies().getMass(part.getID());
  double factor = sqrtT / std::sqrt(mass);

  //Assign the new velocities
  for (size_t iDim = 0; iDim < NDIM; iDim++)
    const_cast<Particle&>(part).getVelocity()[iDim] 
      = Sim->normal_sampler() * factor;

  tmpDat.setDeltaKE(0.5 * mass * (part.getVelocity().nrm2() - tmpDat.getOldVel().nrm2()));
  
  return tmpDat;
}

LNewtonian::LNewtonian(dynamo::SimData* tmp):
  Liouvillean(tmp),
  lastAbsoluteClock(-1),
  lastCollParticle1(0),
  lastCollParticle2(0)  
{}

void
LNewtonian::streamParticle(Particle &particle, const double &dt) const
{
  particle.getPosition() += particle.getVelocity() * dt;

  //The Vector copy is required to make sure that the cached
  //orientation doesn't change during calculation
  if (hasOrientationData())
    orientationData[particle.getID()].orientation 
      = Rodrigues(orientationData[particle.getID()].angularVelocity * dt)
      * Vector(orientationData[particle.getID()].orientation); 
}

double 
LNewtonian::getWallCollision(const Particle &part, 
			     const Vector  &wallLoc, 
			     const Vector  &wallNorm) const
{
  Vector  rij = part.getPosition(),
    vel = part.getVelocity();

  Sim->dynamics.BCs().applyBC(rij, vel);

  double rvdot = (vel | wallNorm);

  rij -= wallLoc;

  if (rvdot < 0)
    return  - ((rij | wallNorm) / rvdot);
  
  return HUGE_VAL;
}


ParticleEventData 
LNewtonian::runWallCollision(const Particle &part, 
			     const Vector  &vNorm,
			     const double& e
			     ) const
{
  updateParticle(part);

  ParticleEventData retVal(part, Sim->dynamics.getSpecies(part), WALL);
  
  const_cast<Particle&>(part).getVelocity()
    -= (1+e) * (vNorm | part.getVelocity()) * vNorm;
  
  retVal.setDeltaKE(0.5 * retVal.getSpecies().getMass(part.getID())
		    * (part.getVelocity().nrm2() 
		       - retVal.getOldVel().nrm2()));
  
  return retVal; 
}

ParticleEventData 
LNewtonian::runAndersenWallCollision(const Particle& part, 
				     const Vector & vNorm,
				     const double& sqrtT
				     ) const
{  
  updateParticle(part);

  if (hasOrientationData())
    M_throw() << "Need to implement thermostating of the rotational degrees"
      " of freedom";

  //This gives a completely new random unit vector with a properly
  //distributed Normal component. See Granular Simulation Book
  ParticleEventData tmpDat(part, Sim->dynamics.getSpecies(part), WALL);
 
  double mass = Sim->dynamics.getSpecies(part).getMass(part.getID());

  for (size_t iDim = 0; iDim < NDIM; iDim++)
    const_cast<Particle&>(part).getVelocity()[iDim] 
      = Sim->normal_sampler() * sqrtT / std::sqrt(mass);
  
  const_cast<Particle&>(part).getVelocity() 
    //This first line adds a component in the direction of the normal
    += vNorm * (sqrtT * sqrt(-2.0*log(1.0-Sim->uniform_sampler()) / mass)
		//This removes the original normal component
		-(part.getVelocity() | vNorm));

  tmpDat.setDeltaKE(0.5 * mass * (part.getVelocity().nrm2() - tmpDat.getOldVel().nrm2()));
  
  return tmpDat; 
}

double
LNewtonian::getSquareCellCollision2(const Particle& part, 
				    const Vector & origin, 
				    const Vector & width) const
{
  Vector  rpos(part.getPosition() - origin);
  Vector  vel(part.getVelocity());
  Sim->dynamics.BCs().applyBC(rpos, vel);
  
#ifdef DYNAMO_DEBUG
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    if ((vel[iDim] == 0) && (std::signbit(vel[iDim])))
      M_throw() << "You have negative zero velocities, don't use them.";
#endif 

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
LNewtonian::getSquareCellCollision3(const Particle& part, 
				    const Vector & origin, 
				    const Vector & width) const
{
  Vector  rpos(part.getPosition() - origin);
  Vector  vel(part.getVelocity());

  Sim->dynamics.BCs().applyBC(rpos, vel);

  int retVal(0);
  double time(HUGE_VAL);
  
#ifdef DYNAMO_DEBUG
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    if ((vel[iDim] == 0) && (std::signbit(vel[iDim])))
      M_throw() << "You have negative zero velocities, dont use them."
		<< "\nPlease think of the neighbour lists.";
#endif

  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    {
      double tmpdt = ((vel[iDim] < 0) 
		      ? -rpos[iDim]/vel[iDim] 
		      : (width[iDim]-rpos[iDim]) / vel[iDim]);

      if (tmpdt < time)
	{
	  time = tmpdt;
	  retVal = (vel[iDim] < 0) ? -(iDim+1) : (iDim+1);
	}
    }

  if (((retVal < 0) && (vel[abs(retVal)-1] > 0))
      || ((retVal > 0) && (vel[abs(retVal)-1] < 0)))
    M_throw() << "Found an error! retVal " << retVal
	      << " vel is " << vel[abs(retVal)-1];

  return retVal;
}

bool 
LNewtonian::DSMCSpheresTest(const Particle& p1, 
			    const Particle& p2, 
			    double& maxprob,
			    const double& factor,
			    CPDData& pdat) const
{
  pdat.vij = p1.getVelocity() - p2.getVelocity();

  //Sim->dynamics.BCs().applyBC(pdat.rij, pdat.vij);
  pdat.rvdot = (pdat.rij | pdat.vij);
  
  if (pdat.rvdot > 0)
    return false; //Positive rvdot

  double prob = factor * (-pdat.rvdot);

  if (prob > maxprob)
    maxprob = prob;

  return prob > Sim->uniform_sampler() * maxprob;
}

PairEventData
LNewtonian::DSMCSpheresRun(const Particle& p1, 
			   const Particle& p2, 
			   const double& e,
			   CPDData& pdat) const
{
  updateParticlePair(p1, p2);  

  PairEventData retVal(p1, p2,
		       Sim->dynamics.getSpecies(p1),
		       Sim->dynamics.getSpecies(p2),
		       CORE);
  
  retVal.rij = pdat.rij;
  retVal.rvdot = pdat.rvdot;

  double p1Mass = retVal.particle1_.getSpecies().getMass(p1.getID());
  double p2Mass = retVal.particle2_.getSpecies().getMass(p2.getID());
  double mu = p1Mass * p2Mass/(p1Mass+p2Mass);

  retVal.dP = retVal.rij * ((1.0 + e) * mu * retVal.rvdot 
			    / retVal.rij.nrm2());  

  //This function must edit particles so it overrides the const!
  const_cast<Particle&>(p1).getVelocity() -= retVal.dP / p1Mass;
  const_cast<Particle&>(p2).getVelocity() += retVal.dP / p2Mass;

  retVal.particle1_.setDeltaKE(0.5 * p1Mass * (p1.getVelocity().nrm2() 
					       - retVal.particle1_.getOldVel().nrm2()));
  
  retVal.particle2_.setDeltaKE(0.5 * p2Mass * (p2.getVelocity().nrm2() 
					       - retVal.particle2_.getOldVel().nrm2()));

  return retVal;
}

PairEventData 
LNewtonian::SmoothSpheresColl(const IntEvent& event, const double& e,
			      const double&, const EEventType& eType) const
{
  const Particle& particle1 = Sim->particleList[event.getParticle1ID()];
  const Particle& particle2 = Sim->particleList[event.getParticle2ID()];

  updateParticlePair(particle1, particle2);  

  PairEventData retVal(particle1, particle2,
		       Sim->dynamics.getSpecies(particle1),
		       Sim->dynamics.getSpecies(particle2),
		       eType);
    
  Sim->dynamics.BCs().applyBC(retVal.rij, retVal.vijold);
  
  double p1Mass = retVal.particle1_.getSpecies().getMass(particle1.getID()); 
  double p2Mass = retVal.particle2_.getSpecies().getMass(particle2.getID());
 
  retVal.rvdot = (retVal.rij | retVal.vijold);

  //Treat special cases if one particle has infinite mass
  if ((p1Mass == 0) && (p2Mass != 0))
    //if (!particle1.testState(Particle::DYNAMIC) && particle2.testState(Particle::DYNAMIC))
    {
      retVal.dP = p2Mass * retVal.rij * ((1.0 + e) * retVal.rvdot / retVal.rij.nrm2());  
      //This function must edit particles so it overrides the const!
      const_cast<Particle&>(particle2).getVelocity() += retVal.dP / p2Mass;
    }
  else 
    if ((p1Mass != 0) && (p2Mass == 0))
      //if (particle1.testState(Particle::DYNAMIC) && !particle2.testState(Particle::DYNAMIC))
      {
	retVal.dP = p1Mass * retVal.rij * ((1.0 + e) * retVal.rvdot / retVal.rij.nrm2());
	//This function must edit particles so it overrides the const!
	const_cast<Particle&>(particle1).getVelocity() -= retVal.dP / p1Mass;
      }
    else
      {
	bool isInfInf = (p1Mass == 0) && (p2Mass == 0);

	//If both particles have infinite mass we just collide them as identical masses
	if (isInfInf) p1Mass = p2Mass = 1;

	double mu = p1Mass * p2Mass / (p1Mass + p2Mass);

	retVal.dP = retVal.rij * ((1.0 + e) * mu * retVal.rvdot / retVal.rij.nrm2());  

	//This function must edit particles so it overrides the const!
	const_cast<Particle&>(particle1).getVelocity() -= retVal.dP / p1Mass;
	const_cast<Particle&>(particle2).getVelocity() += retVal.dP / p2Mass;

	//If both particles have infinite mass we pretend no momentum was transferred
	retVal.dP *= !isInfInf;
      }

  retVal.particle1_.setDeltaKE(0.5 * p1Mass * (particle1.getVelocity().nrm2()
					       - retVal.particle1_.getOldVel().nrm2()));
  
  retVal.particle2_.setDeltaKE(0.5 * p2Mass * (particle2.getVelocity().nrm2() 
					       - retVal.particle2_.getOldVel().nrm2()));

  lastCollParticle1 = particle1.getID();
  lastCollParticle2 = particle2.getID();
  lastAbsoluteClock = Sim->dSysTime;

  return retVal;
}

PairEventData 
LNewtonian::parallelCubeColl(const IntEvent& event, const double& e,
			     const double&, const Matrix& rot,
			     const EEventType& eType) const
{
  const Particle& particle1 = Sim->particleList[event.getParticle1ID()];
  const Particle& particle2 = Sim->particleList[event.getParticle2ID()];

  updateParticlePair(particle1, particle2);

  PairEventData retVal(particle1, particle2,
		       Sim->dynamics.getSpecies(particle1),
		       Sim->dynamics.getSpecies(particle2),
		       eType);
    
  Sim->dynamics.BCs().applyBC(retVal.rij, retVal.vijold);
  
  retVal.rij = rot * Vector(retVal.rij);
  retVal.vijold = rot * Vector(retVal.vijold);

  size_t dim(0);
   
  for (size_t iDim(1); iDim < NDIM; ++iDim)
    if (fabs(retVal.rij[dim]) < fabs(retVal.rij[iDim])) dim = iDim;

  double p1Mass = retVal.particle1_.getSpecies().getMass(particle1.getID()); 
  double p2Mass = retVal.particle2_.getSpecies().getMass(particle2.getID());
  double mu = p1Mass * p2Mass/ (p1Mass + p2Mass);
  
  Vector collvec(0,0,0);

  if (retVal.rij[dim] < 0)
    collvec[dim] = -1;
  else
    collvec[dim] = 1;

  retVal.rvdot = (retVal.rij | retVal.vijold);

  retVal.dP = collvec * (1.0 + e) * mu * (collvec | retVal.vijold);  

  retVal.dP = Transpose(rot) * Vector(retVal.dP);
  retVal.rij = Transpose(rot) * Vector(retVal.rij);
  retVal.vijold = Transpose(rot) * Vector(retVal.vijold);

  //This function must edit particles so it overrides the const!
  const_cast<Particle&>(particle1).getVelocity() -= retVal.dP / p1Mass;
  const_cast<Particle&>(particle2).getVelocity() += retVal.dP / p2Mass;

  retVal.particle1_.setDeltaKE(0.5 * p1Mass * (particle1.getVelocity().nrm2() 
					       - retVal.particle1_.getOldVel().nrm2()));
  
  retVal.particle2_.setDeltaKE(0.5 * p2Mass * (particle2.getVelocity().nrm2() 
					       - retVal.particle2_.getOldVel().nrm2()));

  return retVal;
}

NEventData 
LNewtonian::multibdyCollision(const CRange& range1, const CRange& range2, 
			      const double&, const EEventType& eType) const
{
  Vector COMVel1(0,0,0), COMVel2(0,0,0), COMPos1(0,0,0), COMPos2(0,0,0);
  
  double structmass1(0), structmass2(0);
  
  BOOST_FOREACH(const size_t& ID, range1)
    {
      updateParticle(Sim->particleList[ID]);
      
      double mass = Sim->dynamics.getSpecies(Sim->particleList[ID]).getMass(ID);
      structmass1 += mass;
      
      Vector pos(Sim->particleList[ID].getPosition()),
	vel(Sim->particleList[ID].getVelocity());

      Sim->dynamics.BCs().applyBC(pos, vel);

      COMVel1 += vel * mass;

      COMPos1 += pos * mass;
    }
  
  BOOST_FOREACH(const size_t& ID, range2)
    {
      updateParticle(Sim->particleList[ID]);

      double mass = Sim->dynamics.getSpecies(Sim->particleList[ID]).getMass(ID);
      structmass2 += mass;
      
      Vector pos(Sim->particleList[ID].getPosition()),
	vel(Sim->particleList[ID].getVelocity());

      Sim->dynamics.BCs().applyBC(pos, vel);

      COMVel2 += vel * mass;      

      COMPos2 += pos * mass;
    }
  
  COMVel1 /= structmass1;
  COMVel2 /= structmass2;
  
  COMPos1 /= structmass1;
  COMPos2 /= structmass2;
  
  Vector  rij = COMPos1 - COMPos2, vij = COMVel1 - COMVel2;
  Sim->dynamics.BCs().applyBC(rij, vij);
  double rvdot = (rij | vij);

  double mu = structmass1 * structmass2 / (structmass1 + structmass2);

  static const double e = 1.0;
  Vector  dP = rij * ((1.0 + e) * mu * rvdot / rij.nrm2());

  NEventData retVal;
  BOOST_FOREACH(const size_t& ID, range1)
    {
      ParticleEventData tmpval
	(Sim->particleList[ID],
	 Sim->dynamics.getSpecies(Sim->particleList[ID]),
	 eType);

      const_cast<Particle&>(tmpval.getParticle()).getVelocity()
	-= dP / structmass1;
      
      tmpval.setDeltaKE(0.5 * tmpval.getSpecies().getMass(ID)
			* (tmpval.getParticle().getVelocity().nrm2() 
			   - tmpval.getOldVel().nrm2()));
      
      retVal.L1partChanges.push_back(tmpval);
    }

  BOOST_FOREACH(const size_t& ID, range2)
    {
      ParticleEventData tmpval
	(Sim->particleList[ID],
	 Sim->dynamics.getSpecies(Sim->particleList[ID]),
	 eType);

      const_cast<Particle&>(tmpval.getParticle()).getVelocity()
	+= dP / structmass2;
      
      tmpval.setDeltaKE(0.5 * tmpval.getSpecies().getMass(ID)
			* (tmpval.getParticle().getVelocity().nrm2() 
			   - tmpval.getOldVel().nrm2()));
  
      retVal.L1partChanges.push_back(tmpval);
    }
  
  return retVal;
}

NEventData 
LNewtonian::multibdyWellEvent(const CRange& range1, const CRange& range2, 
			      const double&, const double& deltaKE, EEventType& eType) const
{
  Vector  COMVel1(0,0,0), COMVel2(0,0,0), COMPos1(0,0,0), COMPos2(0,0,0);
  
  double structmass1(0), structmass2(0);
  
  BOOST_FOREACH(const size_t& ID, range1)
    {
      updateParticle(Sim->particleList[ID]);
      double mass = Sim->dynamics.getSpecies(Sim->particleList[ID]).getMass(ID);

      structmass1 += mass;

      Vector pos(Sim->particleList[ID].getPosition()),
	vel(Sim->particleList[ID].getVelocity());
      
      Sim->dynamics.BCs().applyBC(pos, vel);

      COMVel1 += vel * mass;

      COMPos1 += pos * mass;
    }
  
  BOOST_FOREACH(const size_t& ID, range2)
    {
      updateParticle(Sim->particleList[ID]);

      double mass = Sim->dynamics.getSpecies(Sim->particleList[ID]).getMass(ID);
      
      structmass2 += mass;
      
      Vector pos(Sim->particleList[ID].getPosition()),
	vel(Sim->particleList[ID].getVelocity());

      Sim->dynamics.BCs().applyBC(pos, vel);
      COMVel2 += vel * mass;
      COMPos2 += pos * mass;
    }
  
  COMVel1 /= structmass1;
  COMVel2 /= structmass2;
  
  COMPos1 /= structmass1;
  COMPos2 /= structmass2;
  
  Vector  rij = COMPos1 - COMPos2, vij = COMVel1 - COMVel2;
  Sim->dynamics.BCs().applyBC(rij, vij);
  double rvdot = (rij | vij);

  double mu = structmass1 * structmass2 / (structmass1 + structmass2);

  double R2 = rij.nrm2();
  double sqrtArg = rvdot * rvdot + 2.0 * R2 * deltaKE / mu;

  Vector  dP(0,0,0);

  if ((deltaKE < 0) && (sqrtArg < 0))
    {
      eType = BOUNCE;
      dP = rij * 2.0 * mu * rvdot / R2;
    }
  else
    {
      if (deltaKE < 0)
	eType = WELL_KEDOWN;
      else
	eType = WELL_KEUP;
	  
      if (rvdot < 0)
	dP = rij 
	  * (2.0 * deltaKE / (std::sqrt(sqrtArg) - rvdot));
      else
	dP = rij 
	  * (-2.0 * deltaKE / (rvdot + std::sqrt(sqrtArg)));
    }
  
  NEventData retVal;
  BOOST_FOREACH(const size_t& ID, range1)
    {
      ParticleEventData tmpval
	(Sim->particleList[ID],
	 Sim->dynamics.getSpecies(Sim->particleList[ID]),
	 eType);

      const_cast<Particle&>(tmpval.getParticle()).getVelocity()
	-= dP / structmass1;
      
      tmpval.setDeltaKE(0.5 * tmpval.getSpecies().getMass(ID)
			* (tmpval.getParticle().getVelocity().nrm2() 
			   - tmpval.getOldVel().nrm2()));
        
      retVal.L1partChanges.push_back(tmpval);
    }

  BOOST_FOREACH(const size_t& ID, range2)
    {
      ParticleEventData tmpval
	(Sim->particleList[ID],
	 Sim->dynamics.getSpecies(Sim->particleList[ID]),
	 eType);

      const_cast<Particle&>(tmpval.getParticle()).getVelocity()
	+= dP / structmass2;
      
      tmpval.setDeltaKE(0.5 * tmpval.getSpecies().getMass(ID)
			* (tmpval.getParticle().getVelocity().nrm2() 
			   - tmpval.getOldVel().nrm2()));
      
      retVal.L1partChanges.push_back(tmpval);
    }
  
  return retVal;
}

PairEventData 
LNewtonian::SphereWellEvent(const IntEvent& event, const double& deltaKE, 
			    const double &) const
{
  const Particle& particle1 = Sim->particleList[event.getParticle1ID()];
  const Particle& particle2 = Sim->particleList[event.getParticle2ID()];

  updateParticlePair(particle1, particle2);  

  PairEventData retVal(particle1, particle2,
		       Sim->dynamics.getSpecies(particle1),
		       Sim->dynamics.getSpecies(particle2),
		       event.getType());
    
  Sim->dynamics.BCs().applyBC(retVal.rij,retVal.vijold);
  
  retVal.rvdot = (retVal.rij | retVal.vijold);
  
  double p1Mass = retVal.particle1_.getSpecies().getMass(particle1.getID());
  double p2Mass = retVal.particle2_.getSpecies().getMass(particle2.getID());
  double mu = p1Mass * p2Mass / (p1Mass + p2Mass);  
  double R2 = retVal.rij.nrm2();
  double sqrtArg = retVal.rvdot * retVal.rvdot + 2.0 * R2 * deltaKE / mu;
  
  if ((deltaKE < 0) && (sqrtArg < 0))
    {
      event.setType(BOUNCE);
      retVal.setType(BOUNCE);
      retVal.dP = retVal.rij * 2.0 * mu * retVal.rvdot / R2;
    }
  else if (deltaKE==0)
    {
      event.setType(NON_EVENT);
      retVal.setType(NON_EVENT);
      retVal.dP = Vector(0,0,0);
    }
  else
    {
      if (deltaKE < 0)
	{
	  event.setType(WELL_KEDOWN);
	  retVal.setType(WELL_KEDOWN);
	}
      else
	{
	  event.setType(WELL_KEUP);
	  retVal.setType(WELL_KEUP);	  
	}
	  
      retVal.particle1_.setDeltaU(-0.5 * deltaKE);
      retVal.particle2_.setDeltaU(-0.5 * deltaKE);	  
      
      if (retVal.rvdot < 0)
	retVal.dP = retVal.rij 
	  * (2.0 * deltaKE / (std::sqrt(sqrtArg) - retVal.rvdot));
      else
	retVal.dP = retVal.rij 
	  * (-2.0 * deltaKE / (retVal.rvdot + std::sqrt(sqrtArg)));
    }
  
#ifdef DYNAMO_DEBUG
  if (boost::math::isnan(retVal.dP[0]))
    M_throw() << "A nan dp has ocurred";
#endif
  
  //This function must edit particles so it overrides the const!
  const_cast<Particle&>(particle1).getVelocity() -= retVal.dP / p1Mass;
  const_cast<Particle&>(particle2).getVelocity() += retVal.dP / p2Mass;
  
  retVal.particle1_.setDeltaKE(0.5 * p1Mass
			       * (particle1.getVelocity().nrm2() 
				  - retVal.particle1_.getOldVel().nrm2()));
  
  retVal.particle2_.setDeltaKE(0.5 * p2Mass
			       * (particle2.getVelocity().nrm2() 
				  - retVal.particle2_.getOldVel().nrm2()));

  return retVal;
}

void 
LNewtonian::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Type") 
      << "Newtonian";
}

double 
LNewtonian::getPBCSentinelTime(const Particle& part, const double& lMax) const
{
#ifdef DYNAMO_DEBUG
  if (!isUpToDate(part))
    M_throw() << "Particle is not up to date";
#endif

  Vector pos(part.getPosition()), vel(part.getVelocity());

  Sim->dynamics.BCs().applyBC(pos, vel);

  double retval = HUGE_VAL;

  for (size_t i(0); i < NDIM; ++i)
    {
      double tmp = (0.5 * Sim->primaryCellSize[i] - lMax) / fabs(vel[i]);

      if (tmp < retval)
	retval = tmp;
    }

  return retval;
}

std::pair<bool,double>
LNewtonian::getPointPlateCollision(const Particle& part, const Vector& nrw0,
				   const Vector& nhat, const double& Delta,
				   const double& Omega, const double& Sigma,
				   const double& t, bool lastpart) const
{
#ifdef DYNAMO_DEBUG
  if (!isUpToDate(part))
    M_throw() << "Particle1 " << part.getID() << " is not up to date";
#endif
  
  Vector pos(part.getPosition() - nrw0), vel(part.getVelocity());
  Sim->dynamics.BCs().applyBC(pos, vel);

  double t_high;
  double surfaceOffset = pos | nhat;
  double surfaceVel = vel | nhat;
  
  if (surfaceVel > 0)
    t_high = (Sigma + Delta - surfaceOffset) / surfaceVel;
  else
    t_high = -(Sigma + Delta + surfaceOffset) / surfaceVel;
  

  COscillatingPlateFunc fL(vel, nhat, pos, t, Delta, Omega, Sigma);
  
#ifdef DYNAMO_DEBUG
  if (Sigma < 0) M_throw() << "Assuming a positive Sigma here";
#endif

  //A particle has penetrated the plate, probably due to some small numerical error
  //We can just adjust the seperation vector till the particle is on the surface of the plate
  if (fL.F_zeroDeriv() > 0)
    {
#ifdef DYNAMO_DEBUG
      I_cerr() << "Particle is penetrating the \"upper\" plate"
	       << "\nTo avoid rediscovering the root we're adjusting the relative position vector to just touching."
	       << "\nThis is fine if it is a rare event.";
#endif
      fL.fixFZeroSign(false);

#ifdef DYNAMO_DEBUG
      //This is just incase the oscillating plate shape function is broken
      if (fL.F_zeroDeriv() > 0)
	M_throw() << "Failed to adjust the plate position";
#endif
    }
      
  double t_low1 = 0, t_low2 = 0;
  if (lastpart)
    {
      if (-fL.F_zeroDeriv() < fL.F_zeroDerivFlip())
	//Shift the lower bound up so we don't find the same root again
	t_low1 = fabs(2.0 * fL.F_firstDeriv())
	  / fL.F_secondDeriv_max();
      else
	t_low2 = fabs(2.0 * fL.F_firstDeriv())
	  / fL.F_secondDeriv_max();
    }


  //Must be careful with collisions at the end of the interval
  t_high *= 1.01;
  
  std::pair<bool,double> root1 
    = frenkelRootSearch(fL, t_low1, t_high, 1e-12 * Sigma);

  fL.flipSigma();

  if (fL.F_zeroDeriv() < 0)
    {
#ifdef DYNAMO_DEBUG
      I_cerr() << "Particle is penetrating the \"lower\" plate"
	       << "\nTo avoid rediscovering the root we're adjusting the relative position vector to just touching."
	       << "\nThis is fine if it is a rare event.";
#endif
      fL.fixFZeroSign(true);

#ifdef DYNAMO_DEBUG
      //This is just incase the oscillating plate shape function is broken
      if (fL.F_zeroDeriv() < 0)
	M_throw() << "Failed to adjust the plate position";
#endif
    }

  std::pair<bool,double> root2 
    = frenkelRootSearch(fL, t_low2, t_high, 1e-12 * Sigma);

  //Check if the particle is penetrating a wall
  //Or if no roots are found at all
  if ((fabs(surfaceOffset - (nhat | fL.wallPosition())) > Sigma)
      || ((root1.second == HUGE_VAL) && (root2.second == HUGE_VAL))
      || ((t_low1 > t_high) && (t_low2 > t_high)))
    {
      //This can be a problem
#ifdef DYNAMO_DEBUG      
      I_cerr() << "Particle " << part.getID() 
	       << " may be outside/heading out of the plates"
	       << "\nerror = "
	       << (fabs(surfaceOffset - (nhat | fL.wallPosition())) - Sigma) 
	/ Sim->dynamics.units().unitLength()
	       << "\n Root1 = " 
	       << root1.second / Sim->dynamics.units().unitTime()
	       << "\n Root2 = " 
	       << root2.second / Sim->dynamics.units().unitTime();
#endif
      
      //If the particle is going out of bounds, collide now
      if (fL.test_root())
	{
#ifdef DYNAMO_DEBUG
	  {
	    COscillatingPlateFunc ftmp(fL);
	    COscillatingPlateFunc ftmp2(fL);
	    ftmp.flipSigma();
	    
	    double fl01(ftmp.F_zeroDeriv());
	    ftmp.stream(t_low1);
	    double flt_low1(ftmp.F_zeroDeriv());
	    ftmp.stream(t_high - t_low1);
	    double flt_high1(ftmp.F_zeroDeriv());
	    
	    double fl02(ftmp2.F_zeroDeriv());
	    ftmp2.stream(t_low2);
	    double flt_low2(ftmp2.F_zeroDeriv());
	    ftmp2.stream(t_high - t_low2);
	    double flt_high2(ftmp2.F_zeroDeriv());
	    
	    I_cerr() << "****Forcing collision"
		     << "\ndSysTime = " << Sim->dSysTime
		     << "\nlNColl = " << Sim->eventCount
		     << "\nlast part = " << (lastpart ? (std::string("True")) : (std::string("False")))
		     << "\nVel = " << part.getVelocity()[0]
		     << "\nPos = " << part.getPosition()[0]
		     << "\nVwall[0] = " << fL.wallVelocity()[0]
		     << "\nRwall[0] = " << fL.wallPosition()[0]
		     << "\nRwall[0]+Sigma = " << fL.wallPosition()[0] + Sigma
		     << "\nRwall[0]-Sigma = " << fL.wallPosition()[0] - Sigma
		     << "\nSigma + Del = " << Sigma+Delta
		     << "\nGood root = " << fL.test_root()
		     << "\nt_low1 = " << t_low1
		     << "\nt_low2 = " << t_low2
		     << "\nt_high = " << t_high
		     << "\nroot1 = " << root1.second
		     << "\nroot2 = " << root2.second
		     << "\nf1(0) = " << fl01
		     << "\nf1(t_low1) = " << flt_low1
		     << "\nf1(t_high) = " << flt_high1
		     << "\nf2(0)_1 = " << fl02
		     << "\nf2(t_low2) = " << flt_low2
		     << "\nf2(t_high) = " << flt_high2
		     << "\nf'(0) =" << fL.F_firstDeriv()
		     << "\nf''(Max) =" << fL.F_secondDeriv_max()
		     << "\nf(x)=" << (pos | nhat)
		     << "+" << (part.getVelocity() | nhat)
		     << " * x - "
		     << Delta 
		     << " * cos(("
		     << t + Sim->dSysTime << "+ x) * "
		     << Omega << ") - "
		     << Sigma
		     << "; set xrange[0:" << t_high << "]; plot f(x)";
	    ;
	  }
#endif
	  return std::pair<bool, double>(true, 0);
	}
      else
	{
	  //The particle and plate are approaching but might not be
	  //before the overlap is fixed, schedule another test later
	  //on
	  double currRoot = HUGE_VAL;

	  if (root1.first)
	    currRoot = root1.second;

	  if (root2.first && (currRoot > root2.second))
	    currRoot = root2.second;
	  //
	  double tmpt = fabs(surfaceVel - fL.velnHatWall());
	  //This next line sets what the recoil velocity should be
	  //We choose the velocity that gives elastic collisions!
	  tmpt += fL.maxWallVel() * 0.002;
	  tmpt /= fL.F_secondDeriv_max();
	  if (tmpt < currRoot)
	    {
#ifdef DYNAMO_DEBUG
	      I_cout() << "Making a fake collision at " << tmpt << "for particle " << part.getID();
#endif

	      return std::pair<bool,double>(true, tmpt);
	    }
#ifdef DYNAMO_DEBUG
	  else
	    I_cout() << "The current root is lower than the fake one";	    
#endif
	}
    }
  
  return (root1.second < root2.second) ? root1 : root2;
}

ParticleEventData 
LNewtonian::runOscilatingPlate
(const Particle& part, const Vector& rw0, const Vector& nhat, double& delta, 
 const double& omega0, const double& sigma, const double& mass, const double& e, 
 double& t, bool strongPlate) const
{
  std::cout.flush();
  updateParticle(part);

  ParticleEventData retVal(part, Sim->dynamics.getSpecies(part), WALL);

  COscillatingPlateFunc fL(part.getVelocity(), nhat, part.getPosition(),
			   t + Sim->dSysTime, delta, omega0, sigma);

  //Should force the particle to the plate surface

  Vector pos(part.getPosition() - fL.wallPosition()), vel(part.getVelocity());

  Sim->dynamics.BCs().applyBC(pos, vel);
  
  double pmass = retVal.getSpecies().getMass(part.getID());
  double mu = (pmass * mass) / (mass + pmass);

  Vector vwall(fL.wallVelocity());

  //  I_cerr() << "Running event for part " << part.getID() <<
  //	   "\ndSysTime = " << Sim->dSysTime << "\nlNColl = " <<
  //	   Sim->lNColl << "\nVel = " << part.getVelocity()[0] <<
  //	   "\nPos = " << part.getPosition()[0] << "\nVwall[0] = " <<
  //	   fL.wallVelocity()[0] << "\nRwall[0] = " <<
  //	   fL.wallPosition()[0] << "\nRwall[0]+sigma = " <<
  //	   fL.wallPosition()[0] + sigma << "\nRwall[0]-sigma = " <<
  //	   fL.wallPosition()[0] - sigma << "\nsigma + Del = " <<
  //	   sigma+delta << "\nf(0)* = " << fL.F_zeroDeriv() << "\nf'(0)
  //	   =" << fL.F_firstDeriv() << "\nf''(Max) =" <<
  //	   fL.F_secondDeriv_max(0) << "\nf(x)=" << pos[0] << "+" <<
  //	   part.getVelocity()[0] << " * x - " << delta << " * cos(("
  //	   << t << "+ x) * " << omega0 << ") - " << sigma;

  
  //Check the root is valid
  if (!fL.test_root())
    {
      double f0 = fL.F_zeroDeriv(), f1 = fL.F_firstDeriv(),
	f2 = fL.F_secondDeriv_max();
      fL.flipSigma();
      
      I_cerr() <<"Particle " << part.getID()
	       << ", is pulling on the oscillating plate!"
	       << "\nRunning event for part " << part.getID()
	       << "\ndSysTime = " << Sim->dSysTime
	       << "\nlNColl = " << Sim->eventCount
	       << "\nVel = " << part.getVelocity()[0]
	       << "\nPos = " << part.getPosition()[0]
	       << "\nVwall[0] = " << fL.wallVelocity()[0]
	       << "\nRwall[0] = " << fL.wallPosition()[0]
	       << "\nRwall[0]+sigma = " << fL.wallPosition()[0] + sigma
	       << "\nRwall[0]-sigma = " << fL.wallPosition()[0] - sigma
	       << "\nGood root " << fL.test_root()
	       << "\nsigma + Del = " << sigma+delta
	       << "\nf1(0)* = " << fL.F_zeroDeriv()
	       << "\nf1'(0) =" << fL.F_firstDeriv()
	       << "\nf1''(Max) =" << fL.F_secondDeriv_max()
	       << "\nf2(0)* = " << f0
	       << "\nf2'(0) =" << f1
	       << "\nf2''(Max) =" << f2
	       << "\nf(x)=" << (pos | nhat)
	       << "+" << (part.getVelocity() | nhat)
	       << " * x - "
	       << delta
	       << " * cos(("
	       << t + Sim->dSysTime << "+ x) * "
	       << omega0 << ") - "
	       << sigma;
      
      return retVal;
    }

  //static size_t elascount(0);

  double inelas = e;

  double rvdot = ((vel - vwall) | nhat);
  if (fabs(rvdot / fL.maxWallVel()) < 0.002)
    {
      /*
	I_cerr() <<"<<!!!>>Particle " << part.getID() 
	<< " gone elastic!\nratio is " << fabs(rvdot / vwall.nrm())
	<< "\nCount is " << ++elascount;
      */
      inelas = 1.0;
      if (fabs(rvdot / fL.maxWallVel()) < 0.001)
	{
	  if (rvdot < 0)
	    rvdot = -fL.maxWallVel() * 0.01;
	  else
	    rvdot = fL.maxWallVel() * 0.01;
	}
    }

  Vector delP =  nhat * mu * (1.0 + inelas) * rvdot;

  const_cast<Particle&>(part).getVelocity() -=  delP / pmass;

  retVal.setDeltaKE(0.5 * pmass
		    * (part.getVelocity().nrm2() 
		       - retVal.getOldVel().nrm2()));
  
  //Don't progress if you want to not change the plate data
  if (strongPlate) return retVal;

  double numerator = -nhat | ((delP / mass) + vwall);

  double reducedt = Sim->dSysTime 
    - 2.0 * M_PI * int(Sim->dSysTime * omega0 / (2.0*M_PI)) / omega0;
  
  double denominator = omega0 * delta * std::cos(omega0 * (reducedt + t));
  

  double newt = std::atan2(numerator, denominator)/ omega0 
    - Sim->dSysTime;
  
  delta *= std::cos(omega0 * (Sim->dSysTime + t)) 
    / std::cos(omega0 * (Sim->dSysTime + newt));
  
  t = newt;

  t -= 2.0 * M_PI * int(t * omega0 / (2.0*M_PI)) / omega0;

  return retVal; 
}

double 
LNewtonian::getCylinderWallCollision(const Particle& part, 
				     const Vector& wallLoc, 
				     const Vector& wallNorm,
				     const double& radius) const
{
  Vector  rij = part.getPosition() - wallLoc,
    vel = part.getVelocity();

  Sim->dynamics.BCs().applyBC(rij, vel);

  rij -= Vector((rij | wallNorm) * wallNorm);

  vel -= Vector((vel | wallNorm) * wallNorm);

  double B = (vel | rij),
    A = vel.nrm2(),
    C = rij.nrm2() - radius * radius;

  double t = (std::sqrt(B*B - A*C) - B) / A;

  if (boost::math::isnan(t))
    return HUGE_VAL;
  else
    return t;
}

ParticleEventData 
LNewtonian::runCylinderWallCollision(const Particle& part, 
				     const Vector& origin,
				     const Vector& vNorm,
				     const double& e
				     ) const
{
  updateParticle(part);

  ParticleEventData retVal(part, Sim->dynamics.getSpecies(part), WALL);
  
  Vector rij =  origin - part.getPosition();

  Sim->dynamics.BCs().applyBC(rij);

  rij -= Vector((rij | vNorm) * vNorm);

  rij /= rij.nrm();

  const_cast<Particle&>(part).getVelocity()
    -= (1+e) * (rij | part.getVelocity()) * rij;
  
  retVal.setDeltaKE(0.5 * retVal.getSpecies().getMass(part.getID())
		    * (part.getVelocity().nrm2() 
		       - retVal.getOldVel().nrm2()));
  
  return retVal; 
}

ParticleEventData 
LNewtonian::runSphereWallCollision(const Particle& part, 
				   const Vector& origin,
				   const double& e
				   ) const
{
  updateParticle(part);

  ParticleEventData retVal(part, Sim->dynamics.getSpecies(part), WALL);
  
  Vector rij =  origin - part.getPosition();

  Sim->dynamics.BCs().applyBC(rij);

  rij /= rij.nrm();

  const_cast<Particle&>(part).getVelocity()
    -= (1+e) * (rij | part.getVelocity()) * rij;
  
  retVal.setDeltaKE(0.5 * retVal.getSpecies().getMass(part.getID())
		    * (part.getVelocity().nrm2() 
		       - retVal.getOldVel().nrm2()));
  
  return retVal; 
}

bool 
LNewtonian::getLineLineCollision(CPDData& PD, const double& length, 
				 const Particle& p1, const Particle& p2) const
{  
#ifdef DYNAMO_DEBUG
  if (!hasOrientationData())
    M_throw() << "Cannot use this function without orientational data";

  if (!isUpToDate(p1))
    M_throw() << "Particle1 " << p1.getID() << " is not up to date";

  if (!isUpToDate(p2))
    M_throw() << "Particle2 " << p2.getID() << " is not up to date";
#endif

  double t_low = 0.0;
  double t_high = PD.dt;
  
  CLinesFunc fL(PD.rij, PD.vij,
		orientationData[p1.getID()].angularVelocity,
		orientationData[p2.getID()].angularVelocity,
		orientationData[p1.getID()].orientation,
		orientationData[p2.getID()].orientation,
		length);
  
  if (((p1.getID() == lastCollParticle1 && p2.getID() == lastCollParticle2)
       || (p1.getID() == lastCollParticle2 && p2.getID() == lastCollParticle1))
      && Sim->dSysTime == lastAbsoluteClock)
    //Shift the lower bound up so we don't find the same root again
    t_low += fabs(2.0 * fL.F_firstDeriv())
      / fL.F_secondDeriv_max();
  
  //Find window delimited by discs
  std::pair<double,double> dtw = fL.discIntersectionWindow();
  
  if(dtw.first > t_low)
    t_low = dtw.first;
  
  if(dtw.second < t_high)
    t_high = dtw.second;
  
  std::pair<bool,double> root = frenkelRootSearch(fL, t_low, t_high, length * 1e-10);

  if (root.first) 
    { 
      PD.dt = root.second;
      return true; 
    }
  else 
    return false;
}

PairEventData 
LNewtonian::runLineLineCollision(const IntEvent& eevent, const double& elasticity, const double& length) const
{
#ifdef DYNAMO_DEBUG
  if (!hasOrientationData())
    M_throw() << "Cannot use this function without orientational data";
#endif

  const Particle& particle1 = Sim->particleList[eevent.getParticle1ID()];
  const Particle& particle2 = Sim->particleList[eevent.getParticle2ID()];

  updateParticlePair(particle1, particle2);  

  PairEventData retVal(particle1, particle2,
		       Sim->dynamics.getSpecies(particle1),
		       Sim->dynamics.getSpecies(particle2),
		       CORE);
  
  Sim->dynamics.BCs().applyBC(retVal.rij, retVal.vijold);

  retVal.rvdot = (retVal.rij | retVal.vijold);

  double KE1before = getParticleKineticEnergy(particle1);
  double KE2before = getParticleKineticEnergy(particle2);

  CLinesFunc fL(retVal.rij, retVal.vijold,
		orientationData[particle1.getID()].angularVelocity,
		orientationData[particle2.getID()].angularVelocity,
		orientationData[particle1.getID()].orientation,
		orientationData[particle2.getID()].orientation,
		length);

  Vector uPerp = fL.getu1() ^ fL.getu2();

  uPerp /= uPerp.nrm();

  std::pair<double, double> cp = fL.getCollisionPoints();

  // \Delta {\bf v}_{imp}
  Vector  vr = retVal.vijold
    + (cp.first * fL.getw1() ^ fL.getu1()) 
    - (cp.second * fL.getw2() ^ fL.getu2());
  
  double mass = retVal.particle1_.getSpecies().getMass(particle1.getID());
  double inertia = retVal.particle1_.getSpecies().getScalarMomentOfInertia(particle1.getID());

  retVal.dP = uPerp
    * (((vr | uPerp) * (1.0 + elasticity))
       / ((2.0 / mass) + ((cp.first * cp.first + cp.second * cp.second) / inertia)));
  
  const_cast<Particle&>(particle1).getVelocity() -= retVal.dP / mass;
  const_cast<Particle&>(particle2).getVelocity() += retVal.dP / mass;

  orientationData[particle1.getID()].angularVelocity 
    -= (cp.first / inertia) * (fL.getu1() ^ retVal.dP);

  orientationData[particle2.getID()].angularVelocity 
    += (cp.second / inertia) * (fL.getu2() ^ retVal.dP);

  retVal.particle1_.setDeltaKE(getParticleKineticEnergy(particle1) - KE1before);
  retVal.particle2_.setDeltaKE(getParticleKineticEnergy(particle2) - KE2before);

  lastCollParticle1 = particle1.getID();
  lastCollParticle2 = particle2.getID();
  lastAbsoluteClock = Sim->dSysTime;

  return retVal;
}


//Here starts my code for offCenterSpheres
bool 
LNewtonian::getOffCenterSphereOffCenterSphereCollision(CPDData& PD, const double& length,  const double& diameter, 
						       const Particle& p1, const Particle& p2) const

{  
#ifdef DYNAMO_DEBUG
  if (!hasOrientationData())
    M_throw() << "Cannot use this function without orientational data";

  if (!isUpToDate(p1))
    M_throw() << "Particle1 " << p1.getID() << " is not up to date";

  if (!isUpToDate(p2))
    M_throw() << "Particle2 " << p2.getID() << " is not up to date";
#endif

  double t_low = 0.0;
  double t_high = PD.dt;
  double tolerance = 1e-16;
  
  CDumbbellsFunc fL1(PD.rij, PD.vij,
		     orientationData[p1.getID()].angularVelocity,
		     orientationData[p2.getID()].angularVelocity,
		     orientationData[p1.getID()].orientation,
		     orientationData[p2.getID()].orientation,length,diameter);
  
  if (((p1.getID() == lastCollParticle1 && p2.getID() == lastCollParticle2)
       || (p1.getID() == lastCollParticle2 && p2.getID() == lastCollParticle1))
      && Sim->dSysTime == lastAbsoluteClock)
    //Shift the lower bound up so we don't find the same root again
    t_low += fabs(2.0 * fL1.F_firstDeriv())
      / fL1.F_secondDeriv_max();
  
  //I_cout()<<"Sphere intersection between "<<t_low <<" and "<< t_high; 
  std::pair<bool,double> root1 = frenkelRootSearch(fL1, t_low, t_high,length*tolerance);


  CDumbbellsFunc fL2(PD.rij, PD.vij,
		     orientationData[p1.getID()].angularVelocity,
		     orientationData[p2.getID()].angularVelocity,
		     -orientationData[p1.getID()].orientation,
		     orientationData[p2.getID()].orientation,length,diameter);
  
  if (((p1.getID() == lastCollParticle1 && p2.getID() == lastCollParticle2)
       || (p1.getID() == lastCollParticle2 && p2.getID() == lastCollParticle1))
      && Sim->dSysTime == lastAbsoluteClock)
    t_low += fabs(2.0 * fL2.F_firstDeriv())
      / fL2.F_secondDeriv_max();
  
  std::pair<bool,double> root2 = frenkelRootSearch(fL2, t_low, t_high,length*tolerance);

  CDumbbellsFunc fL3(PD.rij, PD.vij,
		     orientationData[p1.getID()].angularVelocity,
		     orientationData[p2.getID()].angularVelocity,
		     orientationData[p1.getID()].orientation,
		     -orientationData[p2.getID()].orientation,length,diameter);
  
  if (((p1.getID() == lastCollParticle1 && p2.getID() == lastCollParticle2)
       || (p1.getID() == lastCollParticle2 && p2.getID() == lastCollParticle1))
      && Sim->dSysTime == lastAbsoluteClock)
    t_low += fabs(2.0 * fL3.F_firstDeriv())
      / fL3.F_secondDeriv_max();
  
  std::pair<bool,double> root3 = frenkelRootSearch(fL3, t_low, t_high,length*tolerance);

  CDumbbellsFunc fL4(PD.rij, PD.vij,
		     orientationData[p1.getID()].angularVelocity,
		     orientationData[p2.getID()].angularVelocity,
		     -orientationData[p1.getID()].orientation,
		     -orientationData[p2.getID()].orientation,length,diameter);
  
  if (((p1.getID() == lastCollParticle1 && p2.getID() == lastCollParticle2)
       || (p1.getID() == lastCollParticle2 && p2.getID() == lastCollParticle1))
      && Sim->dSysTime == lastAbsoluteClock)
    t_low += fabs(2.0 * fL4.F_firstDeriv())
      / fL4.F_secondDeriv_max();
  
  std::pair<bool,double> root4 = frenkelRootSearch(fL4, t_low, t_high,length*tolerance);

  if (root1.first || root2.first || root3.first || root4.first )
    {
      
      double roots[4] = {root1.second, root2.second, root3.second, root4.second};  
      //I_cout()<<roots[0]<< " "<< roots[1]<< " "<<roots[2]<< " "<<roots[3]<< " "; 
      
      PD.dt = *std::min_element(roots,roots+4);
      return true;
    }
  else
    return false;
}

PairEventData 
LNewtonian::runOffCenterSphereOffCenterSphereCollision(const IntEvent& eevent, const double& elasticity, const double& length,const double& diameter) const

{
#ifdef DYNAMO_DEBUG
  if (!hasOrientationData())
    M_throw() << "Cannot use this function without orientational data";
#endif

  const Particle& particle1 = Sim->particleList[eevent.getParticle1ID()];
  const Particle& particle2 = Sim->particleList[eevent.getParticle2ID()];

  updateParticlePair(particle1, particle2);  

  PairEventData retVal(particle1, particle2,
		       Sim->dynamics.getSpecies(particle1),
		       Sim->dynamics.getSpecies(particle2),
		       CORE);
  
  Sim->dynamics.BCs().applyBC(retVal.rij, retVal.vijold);

  retVal.rvdot = (retVal.rij | retVal.vijold);
  I_cout()<<"Two sphere collision\n"; 
  double KE1before = getParticleKineticEnergy(particle1);
  double KE2before = getParticleKineticEnergy(particle2);
  //We need to figure out which two orientations are colliding
  std::pair<int,int> sign;
  double min = HUGE_VAL;
  double norm;
  for (int i=0;i<2;i++)
    {
      for(int j=0;j<2;j++)
	{
	  norm = (retVal.rij +  length * 0.5 * pow(-1,i) * orientationData[particle1.getID()].orientation 
		  -  length * 0.5 * pow(-1,j) * orientationData[particle2.getID()].orientation ).nrm();
	  I_cout()<<"norm " << norm << " dr " << norm - diameter;
	  if(norm < (diameter - 1e-10))
	    {
	      M_throw()<<"`Shit man, we overlap " ;
	    }
	  if(norm < min && fabs(norm - diameter)< 10e-10)
	    {
	      sign.first  = i;
	      sign.second = j;
	      min = norm;
	    }
	}
    }
  I_cout()<<"i " << sign.first << " j " << sign.second; 
  
  //Now the pair sign gives the right two particles
  CDumbbellsFunc fL(retVal.rij, retVal.vijold,
		    orientationData[particle1.getID()].angularVelocity,
		    orientationData[particle2.getID()].angularVelocity,
		    pow(-1,sign.first) * orientationData[particle1.getID()].orientation,
		    -pow(-1,sign.second) * orientationData[particle2.getID()].orientation,length,diameter);

  // Now we have the particles in the moment of the collision
  // then apply collision rules
  Vector u1 = pow(-1,sign.first) * orientationData[particle1.getID()].orientation;
  Vector u2 =  pow(-1,sign.second) * orientationData[particle2.getID()].orientation;

  Vector rhat = retVal.rij + u1 * length / 2 - u2 * length / 2;
  rhat /= rhat.nrm();
  u1 /= u1.nrm();
  u2 /= u2.nrm();

  Vector velContac1 = particle1.getVelocity() 
    + ((orientationData[particle1.getID()].angularVelocity) ^ (((u1 * length) + (rhat * diameter)) / 2));
  Vector velContac2 = particle2.getVelocity() 
    + ((orientationData[particle2.getID()].angularVelocity) ^ (((u2 * length) - (rhat * diameter)) / 2));
  
  Vector velContact = velContac1 - velContac2;
  double mass = retVal.particle1_.getSpecies().getMass(particle1.getID());
  //van Zon's Formulas
  //We need the inertia tensor in the lab frame
  Matrix I1(1.0 / 5.0 * mass * diameter * diameter,0,0,
	    0,1.0 / 5.0 * mass * diameter * diameter + 1.0 / 2.0 * mass * length *length,0,
	    0,0, 1.0 / 5.0 * mass * diameter * diameter + 1.0 / 2.0 * mass * length *length);
  Matrix I2(1.0 / 5.0 * mass * diameter * diameter,0,0,
	    0,1.0 / 5.0 * mass * diameter * diameter + 1.0 / 2.0 * mass * length *length,0,
	    0,0, 1.0 / 5.0 * mass * diameter * diameter + 1.0 / 2.0 * mass * length *length);
  Vector n1 = (u1 * length / 2 + rhat * diameter / 2)^rhat;
  Vector n2 = (u2 * length / 2 - rhat * diameter / 2)^rhat;

  Vector a1 = (rhat - ((rhat|u1) * u1))/((rhat - ((rhat|u1) * u1)).nrm());
  Vector b1 = a1 ^ u1;
  Vector a2 = (rhat - ((rhat|u2) * u2))/((rhat - ((rhat|u2) * u2)).nrm());
  Vector b2 = a2 ^ u2;
  b1 /= b1.nrm();
  b2 /= b2.nrm();

  Vector nI1 = (n1 | u1) * u1 + (n1 | a1) * a1 + (n1 | b1) * b1;  
  Vector nI2 = (n2 | u2) * u2 + (n2 | a2) * a2 + (n2 | b2) * b2;  
  
  double dE1 = nI1 | (Inverse(I1) * nI1);
  double dE2 = nI2 | (Inverse(I2) * nI2);

  double a   = 1/(2 * mass) + (dE1 + dE2)/2;
  double b   = velContact | rhat;
  
  double S = b/a;
  
  Vector vr = retVal.vijold;
  //  retVal.dP = normal * ((vr | normal) * (1.0 + elasticity))/a;
  retVal.dP = rhat * S;
  I_cout()<<"Momentum transfer " << S; 
  I_cout()<<"dv " << vr.nrm(); 
  I_cout()<<"dv at contact " << velContact.nrm();
 
  const_cast<Particle&>(particle1).getVelocity() -= retVal.dP / (2 * mass);
  const_cast<Particle&>(particle2).getVelocity() += retVal.dP / (2 * mass);
 

  //Matrix coordinate transformation
  Matrix W1;
  W1.setRow(0,u1);
  W1.setRow(1,a1);
  W1.setRow(2,b1);
  Matrix W2;
  W2.setRow(0,u2);
  W2.setRow(1,a2);
  W2.setRow(2,b2);


  //I_cout()<<"Orientation1 " << (u1).nrm() << " Orientation2 " << (u2).nrm(); 
  //I_cout()<<"Orientation1 " << (a1).nrm() << " Orientation2 " << (a2).nrm(); 
  //I_cout()<<"Orientation1 " << (b1).nrm() << " Orientation2 " << (b2).nrm(); 


  //Rotational energy before
  I_cout()<<"Energy before " <<  2 * KE1before +
    (I1 * (W1 * orientationData[particle1.getID()].angularVelocity)) *  (W1 * orientationData[particle1.getID()].angularVelocity)/2  + 2 * KE2before +
    (I2 * (W2 * orientationData[particle2.getID()].angularVelocity)) *  (W2 * orientationData[particle2.getID()].angularVelocity)/2 ;
  

  orientationData[particle1.getID()].angularVelocity 
    -= S * ((Inverse(W1) * Inverse(I1) * W1) * n1);
  
  orientationData[particle2.getID()].angularVelocity 
    += S * ((Inverse(W2) * Inverse(I2) * W2) * n2);

  I_cout()<<"Energy after  " << 2 * getParticleKineticEnergy(particle1) +
    (I1 * (W1 * orientationData[particle1.getID()].angularVelocity)) *  (W1 * orientationData[particle1.getID()].angularVelocity)/2  + 2 * getParticleKineticEnergy(particle2) +
    (I2 * (W2 * orientationData[particle2.getID()].angularVelocity)) *  (W2 * orientationData[particle2.getID()].angularVelocity)/2 ;
  Vector velContac1b = particle1.getVelocity() 
    + ((orientationData[particle1.getID()].angularVelocity) ^ (((u1 * length) + (rhat * diameter)) / 2));
  Vector velContac2b = particle2.getVelocity() 
    + ((orientationData[particle2.getID()].angularVelocity) ^ (((u2 * length) - (rhat * diameter)) / 2));
  

  Vector velContactb = velContac1b - velContac2b;
  I_cout()<<"Error in contact velocity " << 1 + (velContact | rhat)/(velContactb | rhat);

  //Done with the collision; keeping track of energy
  retVal.particle1_.setDeltaKE(getParticleKineticEnergy(particle1) - KE1before);
  retVal.particle2_.setDeltaKE(getParticleKineticEnergy(particle2) - KE2before);
  //I_cout()<<"delta energy  " << getParticleKineticEnergy(particle1) - KE1before ;
  //I_cout()<<"delta energy  " <<  getParticleKineticEnergy(particle2) - KE2before;
  lastCollParticle1 = particle1.getID();
  lastCollParticle2 = particle2.getID();
  lastAbsoluteClock = Sim->dSysTime;

  return retVal;
}
//Here ends offCenterSpheres

PairEventData 
LNewtonian::RoughSpheresColl(const IntEvent& event, 
			     const double& e, 
			     const double& et, 
			     const double& d2, 
			     const EEventType& eType
			     ) const
{
#ifdef DYNAMO_DEBUG
  if (!hasOrientationData())
    M_throw() << "Cannot use this function without orientational data";
#endif

  const Particle& particle1 = Sim->particleList[event.getParticle1ID()];
  const Particle& particle2 = Sim->particleList[event.getParticle2ID()];

  updateParticlePair(particle1, particle2);  

  PairEventData retVal(particle1, particle2,
		       Sim->dynamics.getSpecies(particle1),
		       Sim->dynamics.getSpecies(particle2),
		       eType);
    
  Sim->dynamics.BCs().applyBC(retVal.rij, retVal.vijold);
  
  double p1Mass = retVal.particle1_.getSpecies().getMass(particle1.getID()); 
  double p2Mass = retVal.particle2_.getSpecies().getMass(particle2.getID());
  double mu = p1Mass * p2Mass/(p1Mass+p2Mass);
  
  retVal.rvdot = (retVal.rij | retVal.vijold);

  //The normal impulse
  retVal.dP = retVal.rij * ((1.0 + e) * mu * retVal.rvdot / retVal.rij.nrm2());

  Vector eijn = retVal.rij / retVal.rij.nrm();

  //Now the tangential impulse
  Vector gij = retVal.vijold - std::sqrt(d2) * 0.5 
    * ((orientationData[particle1.getID()].angularVelocity
	+ orientationData[particle2.getID()].angularVelocity) ^ eijn);
  
  Vector gijt = (eijn ^ gij) ^ eijn;

  double Jbar = retVal.particle1_.getSpecies().getScalarMomentOfInertia(particle1.getID()) 
    / (p1Mass * d2 * 0.25);
  
  retVal.dP += (Jbar * (1-et) / (2*(Jbar + 1))) * gijt;

  double KE1before = getParticleKineticEnergy(particle1);
  double KE2before = getParticleKineticEnergy(particle2);

  //This function must edit particles so it overrides the const!
  const_cast<Particle&>(particle1).getVelocity() -= retVal.dP / p1Mass;
  const_cast<Particle&>(particle2).getVelocity() += retVal.dP / p2Mass;

  Vector angularVchange = (1-et) / (std::sqrt(d2) * (Jbar+1)) * (eijn ^ gijt);
 
  orientationData[particle1.getID()].angularVelocity
    += angularVchange;
  orientationData[particle2.getID()].angularVelocity 
    += angularVchange;

  retVal.particle1_.setDeltaKE(getParticleKineticEnergy(particle1) - KE1before);
  retVal.particle2_.setDeltaKE(getParticleKineticEnergy(particle2) - KE2before);

  return retVal;
}

ParticleEventData 
LNewtonian::runRoughWallCollision(const Particle& part, 
				  const Vector & vNorm,
				  const double& e,
				  const double& et,
				  const double& r
				  ) const
{
#ifdef DYNAMO_DEBUG
  if (!hasOrientationData())
    M_throw() << "Cannot use this function without orientational data";
#endif

  updateParticle(part);

  ParticleEventData retVal(part, Sim->dynamics.getSpecies(part), WALL);

  double KE1before = getParticleKineticEnergy(part);

  double p1Mass = retVal.getSpecies().getMass(part.getID()); 

  double Jbar = retVal.getSpecies().getScalarMomentOfInertia(part.getID())
    / (p1Mass * r * r);

  Vector gij = part.getVelocity() - r
    * (orientationData[part.getID()].angularVelocity ^ vNorm);
  
  Vector gijt = (vNorm ^ gij) ^ vNorm;
  
  const_cast<Particle&>(part).getVelocity()
    -= (1+e) * (vNorm | part.getVelocity()) * vNorm
    + (Jbar * (1-et) / (Jbar + 1)) * gijt;

  Vector angularVchange = (1-et) / (r * (Jbar+1)) * (vNorm ^ gijt);
 
  orientationData[part.getID()].angularVelocity
    += angularVchange;

  retVal.setDeltaKE(getParticleKineticEnergy(part) - KE1before);
  return retVal; 
}
