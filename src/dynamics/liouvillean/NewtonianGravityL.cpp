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
  particle.getPosition() += particle.getVelocity() * dt;
}

Iflt 
LNewtonianGravity::getWallCollision(const CParticle &part, 
			   const Vector  &wallLoc, 
			   const Vector  &wallNorm) const
{
  D_throw() << "Not implemented yet";    
}

Iflt
LNewtonianGravity::getSquareCellCollision2(const CParticle& part, 
				 const Vector & origin, 
				 const Vector & width) const
{
  D_throw() << "Not implemented yet";    
}

size_t
LNewtonianGravity::getSquareCellCollision3(const CParticle& part, 
				 const Vector & origin, 
				 const Vector & width) const
{
  D_throw() << "Not implemented yet";    
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
