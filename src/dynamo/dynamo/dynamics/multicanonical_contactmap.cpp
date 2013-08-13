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

#include <dynamo/dynamics/multicanonical_contactmap.hpp>
#include <dynamo/interactions/intEvent.hpp>
#include <dynamo/2particleEventData.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/BC/BC.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/species/species.hpp>
#include <dynamo/schedulers/sorters/event.hpp>
#include <dynamo/units/units.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

namespace dynamo {
  DynNewtonianMCCMap::DynNewtonianMCCMap(dynamo::Simulation* tmp, const magnet::xml::Node& XML):
    DynNewtonian(tmp)
  {
    _interaction_name = XML.getAttribute("Interaction");

    if (XML.hasNode("Potential"))
      {
	for (magnet::xml::Node map_node = XML.getNode("Potential").fastGetNode("Map"); 
	     map_node.valid(); ++map_node)
	  {
	    double Wval = map_node.getAttribute("W").as<double>();
	    size_t distance = map_node.getAttribute("Distance").as<size_t>();

	    detail::CaptureMap map;
	    for (magnet::xml::Node entry_node = map_node.fastGetNode("Contact"); entry_node.valid(); ++entry_node)
	      map[detail::CaptureMap::key_type(entry_node.getAttribute("ID1").as<size_t>(), entry_node.getAttribute("ID2").as<size_t>())]
		= entry_node.getAttribute("State").as<size_t>();
	    _W.push_back(std::make_pair(map, WData(distance, Wval)));
	  }
      }
  }

  void 
  DynNewtonianMCCMap::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type")
	<< "NewtonianMCCMap"
	<< magnet::xml::attr("Interaction") << _interaction_name
	<< magnet::xml::tag("Potential");
    
    for (const auto entry : _W)
      {
	XML << magnet::xml::tag("Map")
	    << magnet::xml::attr("W") << entry.second._wval
	    << magnet::xml::attr("Distance") << entry.second._distance
	  ;

	for (const auto& val : entry.first)
	  XML << magnet::xml::tag("Contact")
	      << magnet::xml::attr("ID1") << val.first.first
	      << magnet::xml::attr("ID2") << val.first.second
	      << magnet::xml::attr("State") << val.second
	      << magnet::xml::endtag("Contact");
	XML << magnet::xml::endtag("Map");
      }
    
    XML << magnet::xml::endtag("Potential");
  }

  void DynNewtonianMCCMap::initialise()
  {
    DynNewtonian::initialise();

    if (dynamic_cast<const dynamo::EnsembleNVT*>(Sim->ensemble.get()) == NULL)
      M_throw() << "Multi-canonical simulations require an NVT ensemble";
    
    _interaction = std::dynamic_pointer_cast<ICapture>(Sim->interactions[_interaction_name]);
  }


  NEventData 
  DynNewtonianMCCMap::multibdyWellEvent(const IDRange& range1, const IDRange& range2, 
				  const double&, const double& deltaKE, 
				  EEventType& eType) const
  {
    M_throw() << "Not implemented";
  }

  PairEventData 
  DynNewtonianMCCMap::SphereWellEvent(const IntEvent& event, const double& deltaKE, const double &, size_t newstate) const
  {
    Particle& particle1 = Sim->particles[event.getParticle1ID()];
    Particle& particle2 = Sim->particles[event.getParticle2ID()];

    updateParticlePair(particle1, particle2);  

    PairEventData retVal(particle1, particle2,
			 *Sim->species[particle1],
			 *Sim->species[particle2],
			 event.getType());
    
    Sim->BCs->applyBC(retVal.rij,retVal.vijold);
  
    retVal.rvdot = (retVal.rij | retVal.vijold);
  
    double p1Mass = Sim->species[retVal.particle1_.getSpeciesID()]->getMass(particle1.getID());
    double p2Mass = Sim->species[retVal.particle2_.getSpeciesID()]->getMass(particle2.getID());
    double mu = p1Mass * p2Mass / (p1Mass + p2Mass);  
    double R2 = retVal.rij.nrm2();

    //Calculate the deformed energy change of the system (the one used in the dynamics)
    double MCDeltaKE = deltaKE;

    //If there are entries for the current and possible future energy, then take them into account
    detail::CaptureMap contact_map = *_interaction;
    
    //Add the current bias potential
    MCDeltaKE += W(contact_map) * Sim->ensemble->getEnsembleVals()[2];

    //subtract the possible bias potential in the new state
    contact_map[detail::CaptureMap::key_type(particle1, particle2)] = newstate;
    MCDeltaKE -= W(contact_map) * Sim->ensemble->getEnsembleVals()[2];

    //Test if the deformed energy change allows a capture event to occur
    double sqrtArg = retVal.rvdot * retVal.rvdot + 2.0 * R2 * MCDeltaKE / mu;
    if ((MCDeltaKE < 0) && (sqrtArg < 0))
      {
	event.setType(BOUNCE);
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
    if (boost::math::isnan(retVal.impulse[0]))
      M_throw() << "A nan dp has ocurred";
#endif
  
    //This function must edit particles so it overrides the const!
    particle1.getVelocity() -= retVal.impulse / p1Mass;
    particle2.getVelocity() += retVal.impulse / p2Mass;
  
    return retVal;
  }

  void 
  DynNewtonianMCCMap::replicaExchange(Dynamics& oDynamics)
  {
#ifdef DYNAMO_DEBUG
    if (dynamic_cast<const DynNewtonianMCCMap*>(&oDynamics) == NULL)
      M_throw() << "Trying to swap Dynamicss with different derived types!";
#endif

    DynNewtonianMCCMap& ol(static_cast<DynNewtonianMCCMap&>(oDynamics));
    std::swap(_W, ol._W);
  }

  double 
  DynNewtonianMCCMap::W(const detail::CaptureMap& map) const
  {
    /*Iterate over all tether maps, finding the distance between them
      and looking if the tether applies.*/
    size_t applicable_tethers = 0;
    double accumilated_W = 0;

    for (auto tethermap : _W)
      {
	auto il = tethermap.first.begin();
	auto ir = map.begin();
	
	size_t distance = 0;
	while (il != tethermap.first.end() && ir != map.end())
	  {
	    if (il->first < ir->first)
	      {
		++il;
		++distance;
	      }
	    else if (ir->first < il->first)
	      {
		++ir;
		++distance;
	      }
	    else
	      {
		++il;
		++ir;
	      }
	  }

	while (il != tethermap.first.end())
	  {
	    ++il;
	    ++distance;
	  }

	while (ir != map.end())
	  {
	    ++ir;
	    ++distance;
	  }

	if (distance <= tethermap.second._distance)
	  {
	    ++applicable_tethers;
	    accumilated_W += tethermap.second._wval;
	  }
      }

    return accumilated_W / (applicable_tethers + (applicable_tethers==0));
  }

}
