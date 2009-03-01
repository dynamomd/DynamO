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

CLSLLOD::CLSLLOD(DYNAMO::SimData* tmp):
  CLiouvillean(tmp)
{}

void
CLSLLOD::streamParticle(CParticle& particle, const Iflt& dt) const
{
  particle.getVelocity()[0] += particle.getVelocity()[1] * dt;
}

bool 
CLSLLOD::DSMCSpheresTest(const CParticle& p1, 
			 const CParticle& p2, 
			 Iflt& maxprob,
			 const Iflt& factor,
			 CPDData& pdat) const
{
  pdat.vij = p1.getVelocity() - p2.getVelocity();
  pdat.vij[0] -= pdat.rij[1];
  pdat.rvdot = pdat.rij % pdat.vij;
  
  if (!std::signbit(pdat.rvdot))
    return false; //Positive rvdot
  
  Iflt prob = factor * (-pdat.rvdot);
  
  if (prob > maxprob)
    maxprob = prob;

  return prob > Sim->uniform_sampler() * maxprob;
}

C2ParticleData
CLSLLOD::DSMCSpheresRun(const CParticle& p1, 
			 const CParticle& p2, 
			 const Iflt& e,
			 CPDData& pdat) const
{
  updateParticlePair(p1, p2);  

  C2ParticleData retVal(p1, p2,
			Sim->Dynamics.getSpecies(p1),
			Sim->Dynamics.getSpecies(p2),
			CORE);
  
  retVal.vijold = pdat.vij;

  retVal.rij = pdat.rij;
  retVal.rvdot = pdat.rvdot;

  Iflt p1Mass = retVal.particle1_.getSpecies().getMass(); 
  Iflt p2Mass = retVal.particle2_.getSpecies().getMass();
  Iflt mu = p1Mass * p2Mass/(p1Mass+p2Mass);

  retVal.dP = retVal.rij * ((1.0 + e) * mu * retVal.rvdot 
			    / retVal.rij.square());  

  retVal.calcDeltaKE(mu);

  //This function must edit particles so it overrides the const!
  const_cast<CParticle&>(p1).getVelocity() -= retVal.dP / p1Mass;
  const_cast<CParticle&>(p2).getVelocity() += retVal.dP / p2Mass;

  return retVal;
}

void 
CLSLLOD::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") 
      << "SLLOD";
}

CNParticleData 
CLSLLOD::multibdyCollision(const CRange&, const CRange&, const Iflt&, 
			   const EEventType&) const
{
  D_throw() << "Not Implemented";
}

CNParticleData 
CLSLLOD::multibdyWellEvent(const CRange&, const CRange&, 
			   const Iflt&, const Iflt&, 
			   EEventType&) const
{
  D_throw() << "Not Implemented";
}

bool 
CLSLLOD::SphereSphereInRoot(CPDData& dat, const Iflt& d2) const
{
  D_throw() << "Not Implemented";
}
  
bool 
CLSLLOD::SphereSphereOutRoot(CPDData& dat, const Iflt& d2) const
{
  D_throw() << "Not Implemented";
}

bool 
CLSLLOD::sphereOverlap(const CPDData& dat, const Iflt& d2) const
{
  D_throw() << "Not Implemented";
}

C1ParticleData 
CLSLLOD::randomGaussianEvent(const CParticle& part, const Iflt& sqrtT) const
{
  D_throw() << "Not Implemented";
}

Iflt 
CLSLLOD::getWallCollision(const CParticle &part, 
			   const CVector<> &wallLoc, 
			   const CVector<> &wallNorm) const
{
  D_throw() << "Not Implemented";
}


C1ParticleData 
CLSLLOD::runWallCollision(const CParticle &part, 
			   const CVector<> &vNorm,
			   const Iflt& e
			   ) const
{
  D_throw() << "Not Implemented";
}

C1ParticleData 
CLSLLOD::runAndersenWallCollision(const CParticle& part, 
			 const CVector<>& vNorm,
			 const Iflt& sqrtT
			 ) const
{  
  D_throw() << "Not Implemented";
}

Iflt
CLSLLOD::getSquareCellCollision2(const CParticle& part, 
				 const CVector<>& origin, 
				 const CVector<>& width) const
{
  D_throw() << "Not Implemented";
}

size_t
CLSLLOD::getSquareCellCollision3(const CParticle& part, 
				 const CVector<>& origin, 
				 const CVector<>& width) const
{
  D_throw() << "Not Implemented";
}

C2ParticleData 
CLSLLOD::SmoothSpheresColl(const CIntEvent& event, const Iflt& e,
			   const Iflt&, const EEventType& eType) const
{
  D_throw() << "Not Implemented";
}

C2ParticleData 
CLSLLOD::SphereWellEvent(const CIntEvent& event, const Iflt& deltaKE,
			 const Iflt &) const
{
  D_throw() << "Not Implemented";
}
