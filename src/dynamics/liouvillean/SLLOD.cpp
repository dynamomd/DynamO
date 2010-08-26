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

#include "SLLOD.hpp"
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

LSLLOD::LSLLOD(DYNAMO::SimData* tmp):
  Liouvillean(tmp)
{}

void
LSLLOD::streamParticle(Particle& particle, const Iflt& dt) const
{
  particle.getVelocity()[0] += particle.getVelocity()[1] * dt;
}

bool 
LSLLOD::DSMCSpheresTest(const Particle& p1, 
			 const Particle& p2, 
			 Iflt& maxprob,
			 const Iflt& factor,
			 CPDData& pdat) const
{
  pdat.vij = p1.getVelocity() - p2.getVelocity();
  pdat.vij[0] -= pdat.rij[1];
  pdat.rvdot = (pdat.rij | pdat.vij);
  
  if (!std::signbit(pdat.rvdot))
    return false; //Positive rvdot
  
  Iflt prob = factor * (-pdat.rvdot);
  
  if (prob > maxprob)
    maxprob = prob;

  return prob > Sim->uniform_sampler() * maxprob;
}

PairEventData
LSLLOD::DSMCSpheresRun(const Particle& p1, 
			 const Particle& p2, 
			 const Iflt& e,
			 CPDData& pdat) const
{
  updateParticlePair(p1, p2);  

  PairEventData retVal(p1, p2,
			Sim->dynamics.getSpecies(p1),
			Sim->dynamics.getSpecies(p2),
			CORE);
 
  retVal.vijold = pdat.vij;

  retVal.rij = pdat.rij;
  retVal.rvdot = pdat.rvdot;

  Iflt p1Mass = retVal.particle1_.getSpecies().getMass(); 
  Iflt p2Mass = retVal.particle2_.getSpecies().getMass();
  Iflt mu = p1Mass * p2Mass/(p1Mass+p2Mass);

  retVal.dP = retVal.rij * ((1.0 + e) * mu * retVal.rvdot 
			    / retVal.rij.nrm2());  

  //This function must edit particles so it overrides the const!
  const_cast<Particle&>(p1).getVelocity() -= retVal.dP / p1Mass;
  const_cast<Particle&>(p2).getVelocity() += retVal.dP / p2Mass;

  retVal.particle1_.setDeltaKE(0.5 * retVal.particle1_.getSpecies().getMass()
			       * (p1.getVelocity().nrm2() 
				  - retVal.particle1_.getOldVel().nrm2()));
  
  retVal.particle2_.setDeltaKE(0.5 * retVal.particle2_.getSpecies().getMass()
			       * (p2.getVelocity().nrm2() 
				  - retVal.particle2_.getOldVel().nrm2()));

  return retVal;
}

void 
LSLLOD::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Type") 
      << "SLLOD";
}

NEventData 
LSLLOD::multibdyCollision(const CRange&, const CRange&, const Iflt&, 
			   const EEventType&) const
{
  D_throw() << "Not Implemented";
}

NEventData 
LSLLOD::multibdyWellEvent(const CRange&, const CRange&, 
			   const Iflt&, const Iflt&, 
			   EEventType&) const
{
  D_throw() << "Not Implemented";
}

bool 
LSLLOD::SphereSphereInRoot(CPDData& dat, const Iflt& d2) const
{
  D_throw() << "Not Implemented";
}
  
bool 
LSLLOD::SphereSphereOutRoot(CPDData& dat, const Iflt& d2) const
{
  D_throw() << "Not Implemented";
}

bool 
LSLLOD::sphereOverlap(const CPDData& dat, const Iflt& d2) const
{
  D_throw() << "Not Implemented";
}

ParticleEventData 
LSLLOD::randomGaussianEvent(const Particle& part, const Iflt& sqrtT) const
{
  D_throw() << "Not Implemented";
}

Iflt 
LSLLOD::getWallCollision(const Particle &part, 
			   const Vector  &wallLoc, 
			   const Vector  &wallNorm) const
{
  D_throw() << "Not Implemented";
}


ParticleEventData 
LSLLOD::runWallCollision(const Particle &part, 
			   const Vector  &vNorm,
			   const Iflt& e
			   ) const
{
  D_throw() << "Not Implemented";
}

ParticleEventData 
LSLLOD::runAndersenWallCollision(const Particle& part, 
			 const Vector & vNorm,
			 const Iflt& sqrtT
			 ) const
{  
  D_throw() << "Not Implemented";
}

Iflt
LSLLOD::getSquareCellCollision2(const Particle& part, 
				 const Vector & origin, 
				 const Vector & width) const
{
  D_throw() << "Not Implemented";
}

int
LSLLOD::getSquareCellCollision3(const Particle& part, 
				 const Vector & origin, 
				 const Vector & width) const
{
  D_throw() << "Not Implemented";
}

PairEventData 
LSLLOD::SmoothSpheresColl(const IntEvent& event, const Iflt& e,
			   const Iflt&, const EEventType& eType) const
{
  D_throw() << "Not Implemented";
}

PairEventData 
LSLLOD::SphereWellEvent(const IntEvent& event, const Iflt& deltaKE,
			 const Iflt &) const
{
  D_throw() << "Not Implemented";
}
