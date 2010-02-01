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

#include "NewtonMCL.hpp"
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
#include "../../outputplugins/1partproperty/uenergy.hpp"

LNewtonianMC::LNewtonianMC(DYNAMO::SimData* tmp, const XMLNode& XML):
  LNewtonian(tmp),
  EnergyPotentialStep(1)
{
  if (strcmp(XML.getAttribute("Type"),"NewtonianMC"))
    D_throw() << "Attempting to load NewtonianMC from " 
	      << XML.getAttribute("Type")
	      << " entry";
  try 
    {
      if (XML.isAttributeSet("EnergyStep"))
	EnergyPotentialStep = boost::lexical_cast<Iflt>
	  (XML.getAttribute("EnergyStep"));
      
      EnergyPotentialStep /= Sim->dynamics.units().unitEnergy();

      if (XML.hasChild("PotentialDeformation"))
	{
	  
	  const XMLNode xSubNode(XML.getChildNode("PotentialDeformation"));
	  
	  int xml_iter = 0;
	  
	  unsigned long nEnergies = xSubNode.nChildNode("Entry");
	  
	  //Reserve some space for the entries
	  _MCEnergyPotential.rehash(nEnergies);

	  for (unsigned long i = 0; i < nEnergies; ++i)
	    {
	      XMLNode xBrowseNode = xSubNode.getChildNode("Entry", &xml_iter);

	      Iflt key = 
		boost::lexical_cast<Iflt>(xBrowseNode.getAttribute("Energy")) 
		/ Sim->dynamics.units().unitEnergy();

	      key /= EnergyPotentialStep;

	      Iflt entry = 
		boost::lexical_cast<Iflt>(xBrowseNode.getAttribute("Shift"))
		/ Sim->dynamics.units().unitEnergy();
	      
	      _MCEnergyPotential[int(key + 0.5 - (key < 0))] = entry;
	    }
	}
    }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in LNewtonianMC";
    }
}

void 
LNewtonianMC::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") 
      << "NewtonianMC"
      << xmlw::attr("EnergyStep")
      << EnergyPotentialStep
	  * Sim->dynamics.units().unitEnergy()
      << xmlw::tag("PotentialDeformation");

  typedef std::pair<const Iflt,Iflt> locpair;

  BOOST_FOREACH(const locpair& val, _MCEnergyPotential)
    {
      Iflt key = val.first * EnergyPotentialStep 
	* Sim->dynamics.units().unitEnergy();
      
      Iflt entry = val.second * Sim->dynamics.units().unitEnergy();
	      
      XML << xmlw::tag("Entry")
	  << xmlw::attr("Energy") << key
	  << xmlw::attr("Shift") << entry
	  << xmlw::endtag("Entry");
    }
    
  XML << xmlw::endtag("PotentialDeformation");
}



void LNewtonianMC::initialise()
{
  LNewtonian::initialise();

  if (Sim->getOutputPlugin<OPUEnergy>() == NULL)
    D_throw() << "This liouvillean needs the UEnergy plugin";
}


CNParticleData 
LNewtonianMC::multibdyWellEvent(const CRange& range1, const CRange& range2, 
				const Iflt&, const Iflt& deltaKE, 
				EEventType& eType) const
{
  D_throw() << "Not implemented";
}

C2ParticleData 
LNewtonianMC::SphereWellEvent(const CIntEvent& event, const Iflt& deltaKE, 
			      const Iflt &) const
{
  const CParticle& particle1 = Sim->vParticleList[event.getParticle1ID()];
  const CParticle& particle2 = Sim->vParticleList[event.getParticle2ID()];

  updateParticlePair(particle1, particle2);  

  C2ParticleData retVal(particle1, particle2,
			Sim->dynamics.getSpecies(particle1),
			Sim->dynamics.getSpecies(particle2),
			event.getType());
    
  Sim->dynamics.BCs().applyBC(retVal.rij,retVal.vijold);
  
  retVal.rvdot = (retVal.rij | retVal.vijold);
  
  Iflt p1Mass = retVal.particle1_.getSpecies().getMass();
  Iflt p2Mass = retVal.particle2_.getSpecies().getMass();
  Iflt mu = p1Mass * p2Mass / (p1Mass + p2Mass);  
  Iflt R2 = retVal.rij.nrm2();

  Iflt CurrentE = Sim->getOutputPlugin<OPUEnergy>()->getSimU();

  Iflt Key1FloatVal = CurrentE / EnergyPotentialStep;
  int Key1 = int(Key1FloatVal + 0.5 - (Key1FloatVal < 0));

  Iflt Key2FloatVal = (CurrentE - deltaKE) / EnergyPotentialStep;
  int Key2 = int(Key2FloatVal + 0.5 - (Key2FloatVal < 0));

  Iflt MCDeltaKE = deltaKE;

  boost::unordered_map<int, Iflt>::const_iterator iPtr = _MCEnergyPotential.find(Key1);
  if (iPtr != _MCEnergyPotential.end())
    MCDeltaKE -= iPtr->second;

  iPtr = _MCEnergyPotential.find(Key2);

  if (iPtr != _MCEnergyPotential.end())
    MCDeltaKE -= iPtr->second;

  Iflt sqrtArg = retVal.rvdot * retVal.rvdot + 2.0 * R2 * MCDeltaKE / mu;

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
  if (isnan(retVal.dP[0]))
    D_throw() << "A nan dp has ocurred";
#endif
  
  //This function must edit particles so it overrides the const!
  const_cast<CParticle&>(particle1).getVelocity() -= retVal.dP / p1Mass;
  const_cast<CParticle&>(particle2).getVelocity() += retVal.dP / p2Mass;
  
  retVal.particle1_.setDeltaKE(0.5 * retVal.particle1_.getSpecies().getMass()
			       * (particle1.getVelocity().nrm2() 
				  - retVal.particle1_.getOldVel().nrm2()));
  
  retVal.particle2_.setDeltaKE(0.5 * retVal.particle2_.getSpecies().getMass()
			       * (particle2.getVelocity().nrm2() 
				  - retVal.particle2_.getOldVel().nrm2()));

  return retVal;
}

