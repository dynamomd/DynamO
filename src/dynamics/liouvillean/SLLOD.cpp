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
#include "../../base/is_simdata.hpp"
#include "../species/species.hpp"
#include "../../schedulers/sorters/datastruct.hpp"

LSLLOD::LSLLOD(DYNAMO::SimData* tmp):
  Liouvillean(tmp)
{}

void
LSLLOD::streamParticle(Particle& particle, const double& dt) const
{
  if (particle.getState().testState(ParticleState::DYNAMIC))
    particle.getVelocity()[0] += particle.getVelocity()[1] * dt;
}

bool 
LSLLOD::DSMCSpheresTest(const Particle& p1, 
			 const Particle& p2, 
			 double& maxprob,
			 const double& factor,
			 CPDData& pdat) const
{
  pdat.vij = p1.getVelocity() - p2.getVelocity();
  pdat.vij[0] -= pdat.rij[1];
  pdat.rvdot = (pdat.rij | pdat.vij);
  
  if (!std::signbit(pdat.rvdot))
    return false; //Positive rvdot
  
  double prob = factor * (-pdat.rvdot);
  
  if (prob > maxprob)
    maxprob = prob;

  return prob > Sim->uniform_sampler() * maxprob;
}

PairEventData
LSLLOD::DSMCSpheresRun(const Particle& p1, 
			 const Particle& p2, 
			 const double& e,
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

  double p1Mass = retVal.particle1_.getSpecies().getMass(); 
  double p2Mass = retVal.particle2_.getSpecies().getMass();
  double mu = p1Mass * p2Mass/(p1Mass+p2Mass);

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
LSLLOD::multibdyCollision(const CRange&, const CRange&, const double&, 
			   const EEventType&) const
{
  M_throw() << "Not Implemented";
}

NEventData 
LSLLOD::multibdyWellEvent(const CRange&, const CRange&, 
			   const double&, const double&, 
			   EEventType&) const
{
  M_throw() << "Not Implemented";
}

bool 
LSLLOD::SphereSphereInRoot(CPDData& dat, const double& d2) const
{
  M_throw() << "Not Implemented";
}
  
bool 
LSLLOD::SphereSphereOutRoot(CPDData& dat, const double& d2) const
{
  M_throw() << "Not Implemented";
}

bool 
LSLLOD::sphereOverlap(const CPDData& dat, const double& d2) const
{
  M_throw() << "Not Implemented";
}

ParticleEventData 
LSLLOD::randomGaussianEvent(const Particle& part, const double& sqrtT) const
{
  M_throw() << "Not Implemented";
}

double 
LSLLOD::getWallCollision(const Particle &part, 
			   const Vector  &wallLoc, 
			   const Vector  &wallNorm) const
{
  M_throw() << "Not Implemented";
}


ParticleEventData 
LSLLOD::runWallCollision(const Particle &part, 
			   const Vector  &vNorm,
			   const double& e
			   ) const
{
  M_throw() << "Not Implemented";
}

ParticleEventData 
LSLLOD::runAndersenWallCollision(const Particle& part, 
			 const Vector & vNorm,
			 const double& sqrtT
			 ) const
{  
  M_throw() << "Not Implemented";
}

double
LSLLOD::getSquareCellCollision2(const Particle& part, 
				 const Vector & origin, 
				 const Vector & width) const
{
  M_throw() << "Not Implemented";
}

int
LSLLOD::getSquareCellCollision3(const Particle& part, 
				 const Vector & origin, 
				 const Vector & width) const
{
  M_throw() << "Not Implemented";
}

PairEventData 
LSLLOD::SmoothSpheresColl(const IntEvent& event, const double& e,
			   const double&, const EEventType& eType) const
{
  M_throw() << "Not Implemented";
}

PairEventData 
LSLLOD::SphereWellEvent(const IntEvent& event, const double& deltaKE,
			 const double &) const
{
  M_throw() << "Not Implemented";
}
