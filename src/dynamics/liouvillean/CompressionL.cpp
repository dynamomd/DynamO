/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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
#include "../../base/is_exception.hpp"
#include "../../base/is_simdata.hpp"
#include "../species/species.hpp"

CLCompression::CLCompression(DYNAMO::SimData* tmp, Iflt GR):
  CLNewton(tmp),
  growthRate(GR) {}

bool 
CLCompression::SphereSphereInRoot(CPDData& dat, const Iflt& d2) const
{
  Iflt b = dat.rvdot - d2 
    * (growthRate * growthRate * Sim->dSysTime + growthRate);
  
  if (b < -0.0) 
    {
      Iflt a = dat.v2 - growthRate * growthRate * d2;
      Iflt c = dat.r2 - d2 * (1.0 + growthRate * Sim->dSysTime 
			      * (2.0 + growthRate * Sim->dSysTime));
      Iflt arg = (b * b) - a * c;
      
      if (arg > 0.0) 
	{
	  dat.dt = c / (sqrt(arg) - b);
	  return true;
	}
    }
  
  return false;
}
  
bool 
CLCompression::SphereSphereOutRoot(CPDData& dat, const Iflt& d2) const
{
  Iflt a = dat.v2 - growthRate * growthRate * d2;
  Iflt b = dat.rvdot - d2 * (growthRate * growthRate 
			       * Sim->dSysTime + growthRate);
  Iflt c = d2 * (1.0 + growthRate * Sim->dSysTime 
		 * (2.0 + growthRate * Sim->dSysTime)) - dat.r2;

  Iflt arg = (b * b) + a * c;
  
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
CLCompression::sphereOverlap(const CPDData& dat, const Iflt& d2) const
{
  Iflt currd2 = d2 * (1 + 2.0 * Sim->dSysTime * growthRate 
			+ pow(Sim->dSysTime * growthRate,2));
  
  return ((dat.r2 - currd2) < 0.0);
}

void
CLCompression::streamParticle(CParticle &particle, const Iflt &dt) const
{
  particle.getPosition() +=  particle.getVelocity() * dt;
}

C2ParticleData 
CLCompression::SmoothSpheresColl(const CIntEvent& event, const Iflt& e, const Iflt& d2, const EEventType& eType) const
{
  const CParticle& particle1 = event.getParticle1();
  const CParticle& particle2 = event.getParticle2();

  C2ParticleData retVal(particle1, particle2,
			Sim->Dynamics.getSpecies(particle1),
			Sim->Dynamics.getSpecies(particle2),
			eType);

  Sim->Dynamics.BCs().setPBC(retVal.rij, retVal.vijold);
    
  Iflt p1Mass = retVal.particle1_.getSpecies().getMass(); 
  Iflt p2Mass = retVal.particle2_.getSpecies().getMass(); 
  Iflt mu = p1Mass*p2Mass/(p1Mass+p2Mass);    
  Iflt r2 = retVal.rij.square();
  
  retVal.rvdot = retVal.rij % retVal.vijold;

  retVal.dP = retVal.rij * ((1.0 + e) * mu * (retVal.rvdot - growthRate * sqrt(d2 * r2)) / r2);

  retVal.calcDeltaKE(mu);

  const_cast<CParticle&>(particle1).getVelocity() -= retVal.dP / p1Mass;
  const_cast<CParticle&>(particle2).getVelocity() += retVal.dP / p2Mass;
  
  return retVal;
}

C2ParticleData 
CLCompression::SphereWellEvent(const CIntEvent& event, const Iflt& deltaKE, const Iflt& d2) const
{
  updateParticlePair(event.getParticle1(), event.getParticle2());
  
  C2ParticleData retVal(event.getParticle1(), event.getParticle2(),
			Sim->Dynamics.getSpecies(event.getParticle1()),
			Sim->Dynamics.getSpecies(event.getParticle2()),
			event.getType());
    
  Sim->Dynamics.BCs().setPBC(retVal.rij, retVal.vijold);
    
  Iflt p1Mass = retVal.particle1_.getSpecies().getMass();
  Iflt p2Mass = retVal.particle2_.getSpecies().getMass();
  Iflt mu = p1Mass * p2Mass / (p1Mass + p2Mass);  
  CVector<> urij = retVal.rij.unitVector();

  retVal.rvdot = urij % retVal.vijold;

  Iflt sqrtArg = std::pow(retVal.rvdot - growthRate * sqrt(d2), 2)  + (2.0 * deltaKE / mu);
  if (std::signbit(deltaKE) && std::signbit(sqrtArg))
    {
      event.setType(BOUNCE);
      retVal.setType(BOUNCE);

      retVal.dP = urij * (2.0 * mu * (retVal.rvdot - growthRate * sqrt(d2)));
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
    }

  retVal.rvdot *= retVal.rij.length();
  
  retVal.calcDeltaKE(mu);
  
#ifdef DYNAMO_DEBUG
  if (isnan(retVal.dP[0]))
    D_throw() << "A nan dp has ocurred";
#endif
  
  //This function must edit particles so it overrides the const!
  const_cast<CParticle&>(event.getParticle1()).getVelocity() -= retVal.dP / p1Mass;
  const_cast<CParticle&>(event.getParticle2()).getVelocity() += retVal.dP / p2Mass;
  
  return retVal;
}

void 
CLCompression::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") 
      << "Compression";
}
