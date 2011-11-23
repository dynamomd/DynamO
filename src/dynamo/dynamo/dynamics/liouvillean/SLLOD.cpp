/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
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

#include <dynamo/dynamics/liouvillean/SLLOD.hpp>
#include <dynamo/dynamics/interactions/intEvent.hpp>
#include <dynamo/dynamics/2particleEventData.hpp>
#include <dynamo/dynamics/NparticleEventData.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/dynamics/BC/BC.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <dynamo/dynamics/species/species.hpp>
#include <dynamo/schedulers/sorters/datastruct.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
  LSLLOD::LSLLOD(dynamo::SimData* tmp):
    Liouvillean(tmp)
  {}

  void
  LSLLOD::streamParticle(Particle& particle, const double& dt) const
  {
    if (particle.testState(Particle::DYNAMIC))
      particle.getVelocity()[0] += particle.getVelocity()[1] * dt;
  }

  bool 
  LSLLOD::DSMCSpheresTest(const Particle& p1, 
			  const Particle& p2, 
			  double& maxprob,
			  const double& factor,
			  Vector rij) const
  {
    updateParticlePair(p1, p2);  

    Vector vij = p1.getVelocity() - p2.getVelocity();
    vij[0] -= rij[1];
    double rvdot = (rij | vij);
  
    if (rvdot >= 0)
      return false; //Positive rvdot
  
    double prob = factor * (-rvdot);
  
    if (prob > maxprob)
      maxprob = prob;

    return (prob > Sim->uniform_sampler() * maxprob);
  }

  PairEventData
  LSLLOD::DSMCSpheresRun(const Particle& p1, 
			 const Particle& p2, 
			 const double& e,
			 Vector rij) const
  {
    updateParticlePair(p1, p2);  
    Vector vij = p1.getVelocity() - p2.getVelocity();
    vij[0] -= rij[1];
    double rvdot = (rij | vij);

    PairEventData retVal(p1, p2,
			 Sim->dynamics.getSpecies(p1),
			 Sim->dynamics.getSpecies(p2),
			 CORE);
 
    retVal.vijold = vij;

    retVal.rij = rij;
    retVal.rvdot = rvdot;

    double p1Mass = retVal.particle1_.getSpecies().getMass(p1.getID()); 
    double p2Mass = retVal.particle2_.getSpecies().getMass(p2.getID());
    double mu = p1Mass * p2Mass/(p1Mass+p2Mass);

    retVal.dP = rij * ((1.0 + e) * mu * rvdot / rij.nrm2());  

    //This function must edit particles so it overrides the const!
    const_cast<Particle&>(p1).getVelocity() -= retVal.dP / p1Mass;
    const_cast<Particle&>(p2).getVelocity() += retVal.dP / p2Mass;

    retVal.particle1_.setDeltaKE(0.5 * p1Mass
				 * (p1.getVelocity().nrm2() 
				    - retVal.particle1_.getOldVel().nrm2()));
  
    retVal.particle2_.setDeltaKE(0.5 * p2Mass
				 * (p2.getVelocity().nrm2() 
				    - retVal.particle2_.getOldVel().nrm2()));

    return retVal;
  }

  void 
  LSLLOD::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") 
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

  double 
  LSLLOD::SphereSphereInRoot(const Particle& p1, const Particle& p2, double d) const
  {
    M_throw() << "Not Implemented";
  }
    
  double 
  LSLLOD::SphereSphereOutRoot(const Particle& p1, const Particle& p2, double d) const
  {
    M_throw() << "Not Implemented";
  }

  double 
  LSLLOD::sphereOverlap(const Particle&, const Particle&,
			const double&) const
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
}
