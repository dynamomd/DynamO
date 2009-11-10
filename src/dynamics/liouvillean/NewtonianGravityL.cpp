/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "NewtonianGravityL.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../interactions/intEvent.hpp"
#include "../2particleEventData.hpp"
#include "../NparticleEventData.hpp"
#include "../dynamics.hpp"
#include "../BC/BC.hpp"
#include "../../base/is_exception.hpp"
#include "../../base/is_simdata.hpp"
#include "../species/species.hpp"
#include "../../schedulers/sorters/datastruct.hpp"
#include "shapes/frenkelroot.hpp"
#include "shapes/oscillatingplate.hpp"

LNewtonianGravity::LNewtonianGravity(DYNAMO::SimData* tmp):
  LNewtonian(tmp)
{}

void
LNewtonianGravity::streamParticle(CParticle &particle, const Iflt &dt) const
{
  particle.getPosition() += dt * particle.getVelocity();
  particle.getPosition()[GravityDim] += 0.5 * dt * dt * Gravity;
  particle.getVelocity()[GravityDim] += dt * Gravity;
}

Iflt 
LNewtonianGravity::getWallCollision(const CParticle &part, 
				    const Vector  &wallLoc, 
				    const Vector  &wallNorm) const
{
  Vector  rij = part.getPosition(),
    vel = part.getVelocity();

  Sim->dynamics.BCs().applyBC(rij, vel);

  Iflt adot = wallNorm[GravityDim] * Gravity;
  Iflt vdot = vel | wallNorm;
  Iflt rdot = (rij - wallLoc) | wallNorm;

  Iflt arg = vdot * vdot - 2 * rdot * adot;
  
  if (arg > 0)
    {
      Iflt t = -(vdot + ((vdot<0) ? -1: 1) * std::sqrt(arg));
      Iflt x1 = t / adot;
      Iflt x2 = 2 * rdot / t;

      if (adot > 0)
	//The particle is arcing under the plate
	return (x1 < x2) ? x1 : x2 ;
      else
	//The particle is arcing over the plate
	return (x1 < x2) ? x2 : x1;
    }

  return HUGE_VAL;
}

Iflt
LNewtonianGravity::getSquareCellCollision2(const CParticle& part, 
					   const Vector & origin, 
					   const Vector & width) const
{
  Vector rpos(part.getPosition() - origin);
  Vector vel(part.getVelocity());
  Sim->dynamics.BCs().applyBC(rpos, vel);
  
#ifdef DYNAMO_DEBUG
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    if ((vel[iDim] == 0) && (std::signbit(vel[iDim])))
      D_throw() << "You have negative zero velocities, dont use them."
		<< "\nPlease think of the neighbour lists.";
#endif 

  Iflt retVal = HUGE_VAL;

  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    if (iDim == GravityDim)
      {
	Iflt adot = Gravity;
	Iflt vdot = vel[GravityDim];
	Iflt rdot = rpos[GravityDim];

	Iflt arg = vdot * vdot - 2 * rdot * adot;
	
	if (arg > 0)
	  {
	    Iflt t = -(vdot + ((vdot<0) ? -1: 1) * std::sqrt(arg));
	    Iflt x1 = t / adot;
	    Iflt x2 = 2 * rdot / t;
	    
	    if (Gravity > 0)
	      {
		//The particle is arcing over the boundary
		Iflt tmpdt = (x1 < x2) ? x1 : x2;
		if (tmpdt < retVal)
		  retVal = tmpdt;
	      }
	    else
	      {
		Iflt tmpdt = (x1 < x2) ? x2 : x1;
		if (tmpdt < retVal)
		  retVal = tmpdt;
	      }
	  }
	
	rdot = width[iDim]-rpos[iDim];

	arg = vdot * vdot - 2 * rdot * adot;
	
	if (arg > 0)
	  {
	    Iflt t = -(vdot + ((vdot<0) ? -1: 1) * std::sqrt(arg));
	    Iflt x1 = t / adot;
	    Iflt x2 = 2 * rdot / t;
	    
	    if (Gravity < 0)
	      {
		//The particle is arcing over the boundary
		Iflt tmpdt = (x1 < x2) ? x1 : x2;
		if (tmpdt < retVal)
		  retVal = tmpdt;
	      }
	    else
	      {
		//The particle is arcing under the boundary
		Iflt tmpdt = (x1 < x2) ? x2 : x1;
		if (tmpdt < retVal)
		  retVal = tmpdt;
	      }
	  }
      }
    else
      {
	Iflt tmpdt((vel[iDim] < 0)
		   ? -rpos[iDim]/vel[iDim] 
		   : (width[iDim]-rpos[iDim]) / vel[iDim]);
	
	if (tmpdt < retVal)
	  retVal = tmpdt;
      }
  
  return retVal;
}

int
LNewtonianGravity::getSquareCellCollision3(const CParticle& part, 
					   const Vector & origin, 
					   const Vector & width) const
{
  Vector  rpos(part.getPosition() - origin);
  Vector  vel(part.getVelocity());

  Sim->dynamics.BCs().applyBC(rpos, vel);

  size_t retVal(0);
  Iflt time(HUGE_VAL);
  
#ifdef DYNAMO_DEBUG
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    if ((vel[iDim] == 0) && (std::signbit(vel[iDim])))
      D_throw() << "You have negative zero velocities, dont use them."
		<< "\nPlease think of the neighbour lists.";
#endif

  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    if (iDim == GravityDim)
      {
	Iflt adot = Gravity;
	Iflt vdot = vel[GravityDim];
	Iflt rdot = rpos[GravityDim];
	
	Iflt arg = vdot * vdot - 2 * rdot * adot;
	
	if (arg > 0)
	  {
	    Iflt t = -(vdot + ((vdot<0) ? -1: 1) * std::sqrt(arg));
	    Iflt x1 = t / adot;
	    Iflt x2 = 2 * rdot / t;
	    
	    if (Gravity > 0)
	      {
		//The particle is arcing over the boundary
		Iflt tmpdt = (x1 < x2) ? x1 : x2;
		if (tmpdt < time)
		  {
		    time = tmpdt;
		    retVal = -(iDim+1);
		  }
	      }
	    else
	      {
		Iflt tmpdt = (x1 < x2) ? x2 : x1;
		if (tmpdt < time)
		  {
		    time = tmpdt;
		    retVal = -(iDim+1);
		  }
	      }
	  }
	
	rdot = width[iDim]-rpos[iDim];

	arg = vdot * vdot - 2 * rdot * adot;
	
	if (arg > 0)
	  {
	    Iflt t = -(vdot + ((vdot<0) ? -1: 1) * std::sqrt(arg));
	    Iflt x1 = t / adot;
	    Iflt x2 = 2 * rdot / t;
	    
	    if (Gravity < 0)
	      {
		//The particle is arcing over the boundary
		Iflt tmpdt = (x1 < x2) ? x1 : x2;
		if (tmpdt < time)
		  {
		    time = tmpdt;
		    retVal = iDim;
		  }
	      }
	    else
	      {
		//The particle is arcing under the boundary
		Iflt tmpdt = (x1 < x2) ? x2 : x1;
		if (tmpdt < time)
		  {
		    time = tmpdt;
		    retVal = iDim;
		  }
	      }
	  }
      }  
  else
    {
      Iflt tmpdt = ((vel[iDim] < 0) 
		  ? -rpos[iDim]/vel[iDim] 
		  : (width[iDim]-rpos[iDim]) / vel[iDim]);

      if (tmpdt < time)
	{
	  time = tmpdt;
	  retVal = (vel[iDim] < 0) ? -(iDim+1) : iDim+1;
	}
    }

  return retVal;
}

void 
LNewtonianGravity::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") 
      << "NewtonianGravity";
}

Iflt 
LNewtonianGravity::getPBCSentinelTime(const CParticle& part, const Iflt& lMax) const
{
  D_throw() << "Not implemented yet";  
}

Iflt
LNewtonianGravity::getPointPlateCollision(const CParticle& part, const Vector& nrw0,
				 const Vector& nhat, const Iflt& Delta,
				 const Iflt& Omega, const Iflt& Sigma,
				 const Iflt& t, bool lastpart) const
{
  D_throw() << "Not implemented yet";
}

Iflt 
LNewtonianGravity::getCylinderWallCollision(const CParticle& part, 
				   const Vector& wallLoc, 
				   const Vector& wallNorm,
				   const Iflt& radius) const
{
  Vector  rij = part.getPosition() - wallLoc,
    vel = part.getVelocity();

  Sim->dynamics.BCs().applyBC(rij, vel);

  rij -= Vector((rij | wallNorm) * wallNorm);

  vel -= Vector((vel | wallNorm) * wallNorm);

  Iflt B = (vel | rij),
    A = vel.nrm2(),
    C = rij.nrm2() - radius * radius;

  Iflt t = (std::sqrt(B*B - A*C) - B) / A;

  if (std::isnan(t))
    return HUGE_VAL;
  else
    return t;
}
