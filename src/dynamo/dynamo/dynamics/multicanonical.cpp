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

#include <dynamo/dynamics/multicanonical.hpp>
#include <dynamo/2particleEventData.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/BC/BC.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/species/species.hpp>
#include <dynamo/outputplugins/misc.hpp>
#include <dynamo/units/units.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
  DynNewtonianMC::DynNewtonianMC(dynamo::Simulation* tmp, const magnet::xml::Node& XML):
    DynNewtonian(tmp),
    EnergyPotentialStep(1)
  {
    if (XML.hasNode("PotentialDeformation"))
      {
	EnergyPotentialStep 
	  = XML.getNode("PotentialDeformation").getAttribute("EnergyStep").as<double>()
	  / Sim->units.unitEnergy();

	for (magnet::xml::Node node = XML.getNode("PotentialDeformation").findNode("W"); 
	     node.valid(); ++node)
	  {
	    double energy = node.getAttribute("Energy").as<double>() / Sim->units.unitEnergy();	    
	    double Wval = node.getAttribute("Value").as<double>();
		
	    //Here, the Wval needs to be multiplied by kT to turn it
	    //into an Energy, but the Ensemble is not yet initialised,
	    //we must do this conversion later, when we actually use the W val.
	    _W[lrint(energy / EnergyPotentialStep)] = Wval;
	  }
      }
  }

  void 
  DynNewtonianMC::outputXML(magnet::xml::XmlStream& XML) const
  {
    std::unordered_map<int, double> wout = _W;

    XML << magnet::xml::attr("Type")
	<< "NewtonianMC"
	<< magnet::xml::tag("PotentialDeformation")
	<< magnet::xml::attr("EnergyStep")
	<< EnergyPotentialStep * Sim->units.unitEnergy();

    typedef std::pair<const int, double> locpair;

    for (const locpair& val : wout)
      XML << magnet::xml::tag("W")
	  << magnet::xml::attr("Energy")
	  << val.first * EnergyPotentialStep * Sim->units.unitEnergy()
	  << magnet::xml::attr("Value") << val.second
	  << magnet::xml::endtag("W");
    
    XML << magnet::xml::endtag("PotentialDeformation");
  }



  void DynNewtonianMC::initialise()
  {
    DynNewtonian::initialise();

    if (dynamic_cast<const dynamo::EnsembleNVT*>(Sim->ensemble.get()) == NULL)
      M_throw() << "Multi-canonical simulations require an NVT ensemble";
    
    if (!(Sim->getOutputPlugin<OPMisc>()))
      M_throw() << "Multicanonical dynamics requires the Misc plugin";
  }


  NEventData 
  DynNewtonianMC::multibdyWellEvent(const IDRange& range1, const IDRange& range2, 
				  const double&, const double& deltaKE, 
				  EEventType& eType) const
  {
    M_throw() << "Not implemented";
  }

  PairEventData 
  DynNewtonianMC::SphereWellEvent(Event& event, const double& deltaKE, const double &, size_t) const
  {
    Particle& particle1 = Sim->particles[event._particle1ID];
    Particle& particle2 = Sim->particles[event._particle2ID];

    updateParticlePair(particle1, particle2);  

    PairEventData retVal(particle1, particle2, *Sim->species(particle1), *Sim->species(particle2), event._type);
    
    Sim->BCs->applyBC(retVal.rij,retVal.vijold);
  
    retVal.rvdot = (retVal.rij | retVal.vijold);
  
    double p1Mass = Sim->species[retVal.particle1_.getSpeciesID()]->getMass(particle1.getID());
    double p2Mass = Sim->species[retVal.particle2_.getSpeciesID()]->getMass(particle2.getID());
    double mu = p1Mass * p2Mass / (p1Mass + p2Mass);  
    double R2 = retVal.rij.nrm2();

    double CurrentE = Sim->getOutputPlugin<OPMisc>()->getConfigurationalU();

    //Calculate the deformed energy change of the system (the one used in the dynamics)
    double MCDeltaKE = deltaKE;

    //If there are entries for the current and possible future energy, then take them into account
    MCDeltaKE += W(CurrentE) * Sim->ensemble->getEnsembleVals()[2];
    MCDeltaKE -= W(CurrentE - deltaKE) * Sim->ensemble->getEnsembleVals()[2];

    //Test if the deformed energy change allows a capture event to occur
    double sqrtArg = retVal.rvdot * retVal.rvdot + 2.0 * R2 * MCDeltaKE / mu;
    if ((MCDeltaKE < 0) && (sqrtArg < 0))
      {
	event._type = BOUNCE;
	retVal.setType(BOUNCE);
	retVal.impulse = retVal.rij * 2.0 * mu * retVal.rvdot / R2;
      }
    else
      {
	retVal.particle1_.setDeltaU(-0.5 * deltaKE);
	retVal.particle2_.setDeltaU(-0.5 * deltaKE);	  
      
	if (retVal.rvdot < 0)
	  retVal.impulse = retVal.rij 
	    * (2.0 * MCDeltaKE / (std::sqrt(sqrtArg) - retVal.rvdot));
	else
	  retVal.impulse = retVal.rij 
	    * (-2.0 * MCDeltaKE / (retVal.rvdot + std::sqrt(sqrtArg)));
      }
  
#ifdef DYNAMO_DEBUG
    if (std::isnan(retVal.impulse[0]))
      M_throw() << "A nan dp has ocurred";
#endif
  
    //This function must edit particles so it overrides the const!
    particle1.getVelocity() -= retVal.impulse / p1Mass;
    particle2.getVelocity() += retVal.impulse / p2Mass;
  
    return retVal;
  }

  void DynNewtonianMC::replicaExchange(Dynamics& oDynamics)
  {
#ifdef DYNAMO_DEBUG
    if (dynamic_cast<const DynNewtonianMC*>(&oDynamics) == NULL)
      M_throw() << "Trying to swap Dynamics with different derived types!";
#endif

    DynNewtonianMC& ol(static_cast<DynNewtonianMC&>(oDynamics));

    std::swap(EnergyPotentialStep, ol.EnergyPotentialStep);
    std::swap(_W, ol._W);
  }
}
