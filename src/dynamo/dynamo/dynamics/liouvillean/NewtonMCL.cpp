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

#include <dynamo/dynamics/liouvillean/NewtonMCL.hpp>
#include <dynamo/dynamics/interactions/intEvent.hpp>
#include <dynamo/dynamics/2particleEventData.hpp>
#include <dynamo/dynamics/NparticleEventData.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/dynamics/BC/BC.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <dynamo/dynamics/species/species.hpp>
#include <dynamo/schedulers/sorters/event.hpp>
#include <dynamo/dynamics/liouvillean/shapes/oscillatingplate.hpp>
#include <dynamo/outputplugins/1partproperty/uenergy.hpp>
#include <dynamo/outputplugins/0partproperty/intEnergyHist.hpp>
#include <dynamo/dynamics/units/units.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

namespace dynamo {
  LNewtonianMC::LNewtonianMC(dynamo::SimData* tmp, const magnet::xml::Node& XML):
    LNewtonian(tmp),
    EnergyPotentialStep(1)
  {
    if (strcmp(XML.getAttribute("Type"),"NewtonianMC"))
      M_throw() << "Attempting to load NewtonianMC from " 
		<< XML.getAttribute("Type")
		<< " entry";

    if (dynamic_cast<const dynamo::EnsembleNVT*>(Sim->ensemble.get()) == NULL)
      M_throw() << "Multi-canonical simulations require an NVT ensemble";

    try 
      {     
      
      
	if (XML.hasNode("PotentialDeformation"))
	  {
	    EnergyPotentialStep 
	      = XML.getNode("PotentialDeformation").getAttribute("EnergyStep").as<double>()
	      / Sim->dynamics.units().unitEnergy();

	    for (magnet::xml::Node node = XML.getNode("PotentialDeformation").fastGetNode("W"); 
		 node.valid(); ++node)
	      {
		double energy = node.getAttribute("Energy").as<double>() / Sim->dynamics.units().unitEnergy();	    
		double Wval = node.getAttribute("Value").as<double>();
		
		//Here, the Wval needs to be multiplied by kT to turn it
		//into an Energy, but the Ensemble is not yet initialised,
		//we must do this conversion later, when we actually use the W val.
		_W[lrint(energy / EnergyPotentialStep)] = Wval;
	      }
	  }
      }
    catch (boost::bad_lexical_cast &)
      { M_throw() << "Failed a lexical cast in LNewtonianMC"; }
  }

  void 
  LNewtonianMC::outputXML(magnet::xml::XmlStream& XML) const
  {
    boost::unordered_map<int, double> wout = _W;

    try {
      double pluginBinWidth = Sim->getOutputPlugin<OPIntEnergyHist>()->getBinWidth();

      if (pluginBinWidth != EnergyPotentialStep)
	derr << "WARNING! Multicanonical simulations can only improve the MC potential"
	     << "when the IntEnergyHist bin width (" 
	     << pluginBinWidth * Sim->dynamics.units().unitEnergy()
	     << ") and the MC potential bin widths("
	     << EnergyPotentialStep * Sim->dynamics.units().unitEnergy()
	     << ") match!\nCannot improve potential, preserving old potential."
	  ;
      else
	wout = Sim->getOutputPlugin<OPIntEnergyHist>()->getImprovedW();
    } catch (std::exception&)
      {}

    XML << magnet::xml::attr("Type")
	<< "NewtonianMC"
	<< magnet::xml::tag("PotentialDeformation")
	<< magnet::xml::attr("EnergyStep")
	<< EnergyPotentialStep * Sim->dynamics.units().unitEnergy();

    typedef std::pair<const int, double> locpair;

    BOOST_FOREACH(const locpair& val, wout)
      XML << magnet::xml::tag("W")
	  << magnet::xml::attr("Energy")
	  << val.first * EnergyPotentialStep * Sim->dynamics.units().unitEnergy()
	  << magnet::xml::attr("Value") << val.second
	  << magnet::xml::endtag("W");
    
    XML << magnet::xml::endtag("PotentialDeformation");
  }



  void LNewtonianMC::initialise()
  {
    LNewtonian::initialise();
  }


  NEventData 
  LNewtonianMC::multibdyWellEvent(const Range& range1, const Range& range2, 
				  const double&, const double& deltaKE, 
				  EEventType& eType) const
  {
    M_throw() << "Not implemented";
  }

  PairEventData 
  LNewtonianMC::SphereWellEvent(const IntEvent& event, const double& deltaKE, 
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

    double CurrentE = Sim->getOutputPlugin<OPUEnergy>()->getSimU();

    //Calculate the deformed energy change of the system (the one used in the dynamics)
    double MCDeltaKE = deltaKE;

    //If there are entries for the current and possible future energy, then take them into account
    MCDeltaKE += W(CurrentE) * Sim->ensemble->getEnsembleVals()[2];
    MCDeltaKE -= W(CurrentE - deltaKE) * Sim->ensemble->getEnsembleVals()[2];

    //Test if the deformed energy change allows a capture event to occur
    double sqrtArg = retVal.rvdot * retVal.rvdot + 2.0 * R2 * MCDeltaKE / mu;
    if ((MCDeltaKE < 0) && (sqrtArg < 0))
      {
	event.setType(BOUNCE);
	retVal.setType(BOUNCE);
	retVal.dP = retVal.rij * 2.0 * mu * retVal.rvdot / R2;
      }
    else
      {
	if (MCDeltaKE < 0)
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
	    * (2.0 * MCDeltaKE / (std::sqrt(sqrtArg) - retVal.rvdot));
	else
	  retVal.dP = retVal.rij 
	    * (-2.0 * MCDeltaKE / (retVal.rvdot + std::sqrt(sqrtArg)));
      }
  
#ifdef DYNAMO_DEBUG
    if (boost::math::isnan(retVal.dP[0]))
      M_throw() << "A nan dp has ocurred";
#endif
  
    //This function must edit particles so it overrides the const!
    const_cast<Particle&>(particle1).getVelocity() -= retVal.dP / p1Mass;
    const_cast<Particle&>(particle2).getVelocity() += retVal.dP / p2Mass;
  
    retVal.particle1_.setDeltaKE(0.5 * p1Mass * (particle1.getVelocity().nrm2() 
						 - retVal.particle1_.getOldVel().nrm2()));
  
    retVal.particle2_.setDeltaKE(0.5 * p2Mass * (particle2.getVelocity().nrm2() 
						 - retVal.particle2_.getOldVel().nrm2()));

    return retVal;
  }

  void 
  LNewtonianMC::swapSystem(Liouvillean& oLiouvillean)
  {
#ifdef DYNAMO_DEBUG
    if (dynamic_cast<const LNewtonianMC*>(&oLiouvillean) == NULL)
      M_throw() << "Trying to swap Liouvilleans with different derived types!";
#endif

    LNewtonianMC& ol(static_cast<LNewtonianMC&>(oLiouvillean));

    std::swap(EnergyPotentialStep, ol.EnergyPotentialStep);
    std::swap(_W, ol._W);
  }
}
