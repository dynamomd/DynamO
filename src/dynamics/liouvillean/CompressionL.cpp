/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "CompressionL.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../interactions/intEvent.hpp"
#include "../2particleEventData.hpp"
#include "../dynamics.hpp"
#include "../BC/BC.hpp"
#include "../../base/is_simdata.hpp"
#include "../species/species.hpp"

LCompression::LCompression(DYNAMO::SimData* tmp, double GR):
  LNewtonian(tmp),
  growthRate(GR) {}

bool 
LCompression::SphereSphereInRoot(CPDData& dat, const double& d2) const
{
  double b = dat.rvdot - d2 
    * (growthRate * growthRate * Sim->dSysTime + growthRate);
  
  if (b < -0.0) 
    {
      double a = dat.v2 - growthRate * growthRate * d2;
      double c = dat.r2 - d2 * (1.0 + growthRate * Sim->dSysTime 
			      * (2.0 + growthRate * Sim->dSysTime));
      double arg = (b * b) - a * c;
      
      if (arg > 0.0) 
	{
	  dat.dt = c / (sqrt(arg) - b);
	  return true;
	}
    }
  
  return false;
}
  
bool 
LCompression::SphereSphereOutRoot(CPDData& dat, const double& d2) const
{
  double a = dat.v2 - growthRate * growthRate * d2;
  double b = dat.rvdot - d2 * (growthRate * growthRate 
			       * Sim->dSysTime + growthRate);
  double c = d2 * (1.0 + growthRate * Sim->dSysTime 
		 * (2.0 + growthRate * Sim->dSysTime)) - dat.r2;

  double arg = (b * b) + a * c;
  
  if (arg > 0.0) 
    if (a > 0.0)
      {
	if (std::signbit(b))
	  dat.dt = (sqrt(arg) - b) / a;
	else
	  dat.dt = c / (sqrt(arg) + b);

	return true;
      }
  
  return false;
}

bool 
LCompression::sphereOverlap(const CPDData& dat, const double& d2) const
{
  double currd2 = d2 * (1 + 2.0 * Sim->dSysTime * growthRate 
			+ pow(Sim->dSysTime * growthRate,2));
  
  return ((dat.r2 - currd2) < 0.0);
}

void
LCompression::streamParticle(Particle &particle, const double &dt) const
{
  particle.getPosition() +=  particle.getVelocity() * dt;
}

PairEventData 
LCompression::SmoothSpheresColl(const IntEvent& event, const double& e, const double& d2, const EEventType& eType) const
{
  const Particle& particle1 = Sim->particleList[event.getParticle1ID()];
  const Particle& particle2 = Sim->particleList[event.getParticle2ID()];

  updateParticlePair(particle1, particle2);  

  PairEventData retVal(particle1, particle2,
			Sim->dynamics.getSpecies(particle1),
			Sim->dynamics.getSpecies(particle2),
			eType);

  Sim->dynamics.BCs().applyBC(retVal.rij, retVal.vijold);
    
  double p1Mass = retVal.particle1_.getSpecies().getMass(); 
  double p2Mass = retVal.particle2_.getSpecies().getMass(); 
  double mu = p1Mass*p2Mass/(p1Mass+p2Mass);    
  double r2 = retVal.rij.nrm2();
  
  retVal.rvdot = (retVal.rij | retVal.vijold);

  retVal.dP = retVal.rij * ((1.0 + e) * mu * (retVal.rvdot - growthRate * sqrt(d2 * r2)) / r2);

  const_cast<Particle&>(particle1).getVelocity() -= retVal.dP / p1Mass;
  const_cast<Particle&>(particle2).getVelocity() += retVal.dP / p2Mass;

  retVal.particle1_.setDeltaKE(0.5 * retVal.particle1_.getSpecies().getMass()
			       * (particle1.getVelocity().nrm2() 
				  - retVal.particle1_.getOldVel().nrm2()));
  
  retVal.particle2_.setDeltaKE(0.5 * retVal.particle2_.getSpecies().getMass()
			       * (particle2.getVelocity().nrm2()
				  - retVal.particle2_.getOldVel().nrm2()));
  
  return retVal;
}

PairEventData 
LCompression::SphereWellEvent(const IntEvent& event, const double& deltaKE, const double& d2) const
{
  const Particle& particle1 = Sim->particleList[event.getParticle1ID()];
  const Particle& particle2 = Sim->particleList[event.getParticle2ID()];

  updateParticlePair(particle1, particle2);  
  
  PairEventData retVal(particle1, particle2,
			Sim->dynamics.getSpecies(particle1),
			Sim->dynamics.getSpecies(particle2),
			event.getType());
    
  Sim->dynamics.BCs().applyBC(retVal.rij, retVal.vijold);
    
  double p1Mass = retVal.particle1_.getSpecies().getMass();
  double p2Mass = retVal.particle2_.getSpecies().getMass();
  double mu = p1Mass * p2Mass / (p1Mass + p2Mass);  
  Vector  urij = retVal.rij / retVal.rij.nrm();

  retVal.rvdot = (urij | retVal.vijold);

  double sqrtArg = std::pow(retVal.rvdot - growthRate * sqrt(d2), 2)  + (2.0 * deltaKE / mu);
  if (std::signbit(deltaKE) && std::signbit(sqrtArg))
    {
      event.setType(BOUNCE);
      retVal.setType(BOUNCE);

      retVal.dP = urij * (2.0 * mu * (retVal.rvdot - growthRate * sqrt(d2)));
    }
  else if (deltaKE==0)
    {
      event.setType(NON_EVENT);
      retVal.setType(NON_EVENT);
      retVal.dP = Vector(0,0,0);
    }
  else
    {
      if (std::signbit(deltaKE))
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
      
      if (std::signbit(retVal.rvdot))
	retVal.dP = urij 
	  * (2.0 * deltaKE / (growthRate * sqrt(d2) + std::sqrt(sqrtArg) - retVal.rvdot ));
      else
	retVal.dP = urij 
	  * (2.0 * deltaKE / (growthRate * sqrt(d2) - std::sqrt(sqrtArg) - retVal.rvdot ));
	;
    }

  retVal.rvdot *= retVal.rij.nrm();
  
#ifdef DYNAMO_DEBUG
  if (isnan(retVal.dP[0]))
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
  const_cast<Particle&>(particle1).getVelocity() -= retVal.dP / p1Mass;
  const_cast<Particle&>(particle2).getVelocity() += retVal.dP / p2Mass;
  
  retVal.particle1_.setDeltaKE(0.5 * retVal.particle1_.getSpecies().getMass()
			       * (particle1.getVelocity().nrm2() 
				  - retVal.particle1_.getOldVel().nrm2()));
  
  retVal.particle2_.setDeltaKE(0.5 * retVal.particle2_.getSpecies().getMass()
			       * (particle2.getVelocity().nrm2() 
				  - retVal.particle2_.getOldVel().nrm2()));
  
  return retVal;
}

void 
LCompression::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Type") 
      << "Compression";
}

double 
LCompression::getPBCSentinelTime(const Particle& part,
				  const double& lMax) const
{
#ifdef DYNAMO_DEBUG
  if (!isUpToDate(part))
    M_throw() << "Particle is not up to date";
#endif

  Vector  pos(part.getPosition()), vel(part.getVelocity());

  Sim->dynamics.BCs().applyBC(pos, vel);

  double retval = (0.5 * Sim->aspectRatio[0] - lMax) / (fabs(vel[0]) + lMax * growthRate);

  for (size_t i(1); i < NDIM; ++i)
    {
      double tmp = (0.5 * Sim->aspectRatio[i] - lMax) / (fabs(vel[0]) + lMax * growthRate);
      
      if (tmp < retval)
	retval = tmp;
    }

  return retval;
}

PairEventData 
LCompression::parallelCubeColl(const IntEvent&, const double&,
				const double&, const EEventType&) const
{ M_throw() << "Not Implemented"; }
