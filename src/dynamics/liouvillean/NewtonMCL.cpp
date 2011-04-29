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

#include "NewtonMCL.hpp"
#include "../interactions/intEvent.hpp"
#include "../2particleEventData.hpp"
#include "../NparticleEventData.hpp"
#include "../dynamics.hpp"
#include "../BC/BC.hpp"
#include "../../base/is_simdata.hpp"
#include "../species/species.hpp"
#include "../../schedulers/sorters/datastruct.hpp"
#include "shapes/frenkelroot.hpp"
#include "shapes/oscillatingplate.hpp"
#include "../../outputplugins/1partproperty/uenergy.hpp"
#include "../units/units.hpp"
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

LNewtonianMC::LNewtonianMC(dynamo::SimData* tmp, const magnet::xml::Node& XML):
  LNewtonian(tmp),
  EnergyPotentialStep(1)
{
  if (strcmp(XML.getAttribute("Type"),"NewtonianMC"))
    M_throw() << "Attempting to load NewtonianMC from " 
	      << XML.getAttribute("Type")
	      << " entry";
  try 
    {
      if (XML.getAttribute("EnergyStep").valid())
	EnergyPotentialStep = XML.getAttribute("EnergyStep").as<double>();
      
      EnergyPotentialStep /= Sim->dynamics.units().unitEnergy();
      
      if (XML.getNode("PotentialDeformation").valid())
	for (magnet::xml::Node node = XML.getNode("PotentialDeformation").getNode("Entry"); 
	     node.valid(); ++node)
	  {
	    double key = node.getAttribute("Energy").as<double>()
	      / Sim->dynamics.units().unitEnergy();
	    
	    key /= EnergyPotentialStep;
	    _MCEnergyPotential[int(key + 0.5 - (key < 0))] 
	      = node.getAttribute("Shift").as<double>()
	      / Sim->dynamics.units().unitEnergy();
	  }
    }
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in LNewtonianMC";
    }
}

void 
LNewtonianMC::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Type") 
      << "NewtonianMC"
      << xml::attr("EnergyStep")
      << EnergyPotentialStep
	  * Sim->dynamics.units().unitEnergy()
      << xml::tag("PotentialDeformation");

  typedef std::pair<const double,double> locpair;

  BOOST_FOREACH(const locpair& val, _MCEnergyPotential)
    {
      double key = val.first * EnergyPotentialStep 
	* Sim->dynamics.units().unitEnergy();
      
      double entry = val.second * Sim->dynamics.units().unitEnergy();
	      
      XML << xml::tag("Entry")
	  << xml::attr("Energy") << key
	  << xml::attr("Shift") << entry
	  << xml::endtag("Entry");
    }
    
  XML << xml::endtag("PotentialDeformation");
}



void LNewtonianMC::initialise()
{
  LNewtonian::initialise();

  if (Sim->getOutputPlugin<OPUEnergy>() == NULL)
    M_throw() << "This liouvillean needs the UEnergy plugin";
}


NEventData 
LNewtonianMC::multibdyWellEvent(const CRange& range1, const CRange& range2, 
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

  double Key1FloatVal = CurrentE / EnergyPotentialStep;
  int Key1 = int(Key1FloatVal + 0.5 - (Key1FloatVal < 0));

  double Key2FloatVal = (CurrentE - deltaKE) / EnergyPotentialStep;
  int Key2 = int(Key2FloatVal + 0.5 - (Key2FloatVal < 0));

  double MCDeltaKE = deltaKE;

  boost::unordered_map<int, double>::const_iterator iPtr = _MCEnergyPotential.find(Key1);
  if (iPtr != _MCEnergyPotential.end())
    MCDeltaKE -= iPtr->second;

  iPtr = _MCEnergyPotential.find(Key2);

  if (iPtr != _MCEnergyPotential.end())
    MCDeltaKE -= iPtr->second;

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
  
#ifdef dynamo_DEBUG
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

