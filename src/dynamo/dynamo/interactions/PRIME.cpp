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

#include <dynamo/interactions/PRIME.hpp>
#include <dynamo/topology/PRIME.hpp>
#include <dynamo/BC/BC.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/globals/global.hpp>
#include <dynamo/particle.hpp>
#include <dynamo/interactions/intEvent.hpp>
#include <dynamo/species/species.hpp>
#include <dynamo/2particleEventData.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/ranges/IDRangeRange.hpp>
#include <dynamo/ranges/IDPairRangeSingle.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <cmath>
#include <iomanip>
#include <algorithm>

namespace dynamo {
  IPRIME::IPRIME(const magnet::xml::Node& XML, dynamo::Simulation* tmp):
    ICapture(tmp, NULL) //A temporary value!
  {
    operator<<(XML);
  }

  void 
  IPRIME::operator<<(const magnet::xml::Node& XML)
  { 
    Interaction::operator<<(XML);

    const std::string topologyName = XML.getAttribute("Topology");
    _PRIME_HB_strength = XML.getAttribute("HBStrength").as<double>();
    _topology = std::dynamic_pointer_cast<TPRIME>(Sim->topology[topologyName]);
    
    if (!_topology)
      M_throw() << "For \"" << getName() << "\", Topology is not a PRIME topology.";

    ICapture::loadCaptureMap(XML);

    if (XML.hasNode("HBonds"))
      for (magnet::xml::Node node = XML.getNode("HBonds").findNode("Bond"); node.valid(); ++node)
	_HBonds.insert(HbondMapType::value_type(node.getAttribute("NH").as<size_t>(), node.getAttribute("CO").as<size_t>()));
  }

  void 
  IPRIME::initialise(size_t nID)
  {
    Interaction::initialise(nID);
    ICapture::initCaptureMap();

    //Need to initialise the HBond map!
  }

  size_t
  IPRIME::captureTest(const Particle& p1, const Particle& p2) const
  {
    if (&(*(Sim->getInteraction(p1, p2))) != this) return false;

    //Calculate the interaction parameters (and name them sensibly)
    const auto interaction_data = getInteractionParameters(p1.getID(), p2.getID());    
    const double outer_diameter = std::get<0>(interaction_data),
      bond_energy = std::get<2>(interaction_data);
    
    //Check that this has a finite well energy
    if (!std::isfinite(bond_energy)) return false;

#ifdef DYNAMO_DEBUG
    const double inner_diameter = std::get<1>(interaction_data);

    if (Sim->dynamics->sphereOverlap(p1, p2, inner_diameter))
      derr << "Warning! Two particles might be overlapping"
	   << "Overlap is " << Sim->dynamics->sphereOverlap(p1, p2, inner_diameter) / Sim->units.unitLength()
	   << "\nd = " << inner_diameter / Sim->units.unitLength() << std::endl;
#endif
 
    return Sim->dynamics->sphereOverlap(p1, p2, outer_diameter) > 0;
  }

  double 
  IPRIME::getInternalEnergy() const
  { 
    //Once the capture maps are loaded just iterate through that determining energies
    double Energy = 0.0;
    for (const ICapture::value_type& IDs : *this)
      Energy += getInternalEnergy(Sim->particles[IDs.first.first], Sim->particles[IDs.first.second]);
    return Energy; 
  }

  double 
  IPRIME::getInternalEnergy(const Particle& p1, const Particle& p2) const
  {
    const auto interaction_data = getInteractionParameters(p1.getID(), p2.getID());    
    const double bond_energy = std::get<2>(interaction_data);
    return bond_energy * isCaptured(p1, p2);
  }


  std::array<double, 4> 
  IPRIME::getGlyphSize(size_t ID) const
  {
    return {{TPRIME::_PRIME_diameters[getBeadData(ID).bead_type], 0, 0, 0}};
  }

  double
  IPRIME::getExcludedVolume(size_t ID) const
  {
    //This calculation only includes the volumes which are always
    //excluded (i.e. the hard core)
    double diam = TPRIME::_PRIME_diameters[getBeadData(ID).bead_type];
    return diam * diam * diam * M_PI / 6.0;
  }

  double
  IPRIME::maxIntDist() const
  {
    double maxdiam = 0;
    maxdiam = std::max(maxdiam, *std::max_element(TPRIME::_PRIME_diameters, TPRIME::_PRIME_diameters + 3));
    maxdiam = std::max(maxdiam, (1.0 + TPRIME::_PRIME_bond_tolerance) * (*std::max_element(TPRIME::_PRIME_BB_bond_lengths, TPRIME::_PRIME_BB_bond_lengths + 3 * 3)));
    maxdiam = std::max(maxdiam, (1.0 + TPRIME::_PRIME_bond_tolerance) * (*std::max_element(TPRIME::_PRIME_pseudobond_lengths, TPRIME::_PRIME_pseudobond_lengths + 3 * 3)));
    maxdiam = std::max(maxdiam, (1.0 + TPRIME::_PRIME_bond_tolerance) * TPRIME::_PRIME_CH_CH_pseudobond_length);
    return maxdiam;
  }

  std::tuple<double, double, double, size_t, size_t>
  IPRIME::getInteractionParameters(size_t pID1, size_t pID2) const
  {
    TPRIME::BeadData p1Data = getBeadData(pID1);
    TPRIME::BeadData p2Data = getBeadData(pID2);
    
    //Ensure that the first bead has the lowest bead_type (it simplifies the logic later)
    if (p1Data.bead_type > p2Data.bead_type) {
      std::swap(pID1, pID2);
      std::swap(p1Data, p2Data);
    }

    const size_t no_HB_res = std::numeric_limits<size_t>::max();

    if (p1Data.bead_type > TPRIME::CO) //SC-SC interaction (as p2Data.bead_type >= p1Data.bead_type)
      {
        const double inner_diameter = TPRIME::_PRIME_diameters[ 22 * p1Data.bead_type + p2Data.bead_type];
        const double outer_diameter = TPRIME::_PRIME_well_diameters[ 22 * p1Data.bead_type + p2Data.bead_type];
        const double bond_energy    = TPRIME::_PRIME_well_depths[ 22 * p1Data.bead_type + p2Data.bead_type];
	return std::make_tuple(outer_diameter, inner_diameter, bond_energy, no_HB_res, no_HB_res);
      }
    else if (p2Data.bead_type <= TPRIME::CO) //BB-BB interaction (as p1Data.bead_type <= p2Data.bead_type)
      {
        const size_t loc1 = p1Data.bead_type + 3 * p1Data.residue;
        const size_t loc2 = p2Data.bead_type + 3 * p2Data.residue;
        const size_t distance = std::max(loc1, loc2) - std::min(loc1, loc2);

        //This treats the special cases if they are 0,1,2, or three backbone
        //bonds apart
        switch (distance)
          {
          case 0:
            M_throw() << "Invalid backbone distance of 0";
          case 1:
            {
	      //Every type of this interaction is a bonded interaction
              const double inner_diameter = TPRIME::_PRIME_BB_bond_lengths[3 * p1Data.bead_type + p2Data.bead_type] * (1.0 - TPRIME::_PRIME_bond_tolerance);
              const double outer_diameter = TPRIME::_PRIME_BB_bond_lengths[3 * p1Data.bead_type + p2Data.bead_type] * (1.0 + TPRIME::_PRIME_bond_tolerance);
              const double bond_energy    = -std::numeric_limits<double>::infinity();
	      return std::make_tuple(outer_diameter, inner_diameter, bond_energy, no_HB_res, no_HB_res);
            }
            break;
          case 2:
            {
	      //Every type of this interaction is a pseudobond interaction
              const double inner_diameter = TPRIME::_PRIME_pseudobond_lengths[3 * p1Data.bead_type + p2Data.bead_type] * (1.0 - TPRIME::_PRIME_bond_tolerance);
              const double outer_diameter = TPRIME::_PRIME_pseudobond_lengths[3 * p1Data.bead_type + p2Data.bead_type] * (1.0 + TPRIME::_PRIME_bond_tolerance);
              const double bond_energy = -std::numeric_limits<double>::infinity();
	      return std::make_tuple(outer_diameter, inner_diameter, bond_energy, no_HB_res, no_HB_res);
            }
            break;
          case 3:
            {
              //Check if this is the special pseudobond
              if ((p1Data.bead_type == TPRIME::CH) && (p2Data.bead_type == TPRIME::CH))
                {
                  const double inner_diameter = TPRIME::_PRIME_CH_CH_pseudobond_length * (1.0 - TPRIME::_PRIME_bond_tolerance);
                  const double outer_diameter = TPRIME::_PRIME_CH_CH_pseudobond_length * (1.0 + TPRIME::_PRIME_bond_tolerance);
                  const double bond_energy = -std::numeric_limits<double>::infinity();
		  return std::make_tuple(outer_diameter, inner_diameter, bond_energy, no_HB_res, no_HB_res);
                }
              else
		{
		  //Close backbone-backbone hard-sphere interaction
		  const double inner_diameter = 0.0;
		  const double outer_diameter = TPRIME::_PRIME_diameters[22 * p1Data.bead_type + p2Data.bead_type] * TPRIME::_PRIME_near_diameter_scale_factor;
		  const double bond_energy = std::numeric_limits<double>::infinity();
		  return std::make_tuple(outer_diameter, inner_diameter, bond_energy, no_HB_res, no_HB_res);
		}
	    }
            break;
          default:
	    {
	      //Determine if HB interactions are present or not
	      //If the time-independent criteria are met, it's not a hard-sphere and we track distances.
	      if ((p1Data.bead_type == TPRIME::NH) && (p2Data.bead_type == TPRIME::CO) && (p1Data.location != TPRIME::NH_END) && (p2Data.location != TPRIME::CO_END) && (std::abs(int(p1Data.residue) - int(p2Data.residue)) > 3))
		{
		  const size_t NH_res = p1Data.residue;
		  const size_t CO_res = p2Data.residue;
		  const double inner_diameter = TPRIME::_PRIME_diameters[ 22 * p1Data.bead_type + p2Data.bead_type ];
		  const double outer_diameter = TPRIME::_PRIME_HB_well_diameter;
		  const double bond_energy = checkTimeDependentCriteria(NH_res, CO_res, 0) ? -_PRIME_HB_strength : 0;
		  return std::make_tuple(outer_diameter, inner_diameter, bond_energy, NH_res, CO_res);
		}
	      else if ((p1Data.bead_type == TPRIME::CH) && (p2Data.bead_type == TPRIME::CO) && (p1Data.location != TPRIME::NH_END) && (p2Data.location != TPRIME::CO_END) && (std::abs(int(p1Data.residue) - int(p2Data.residue)) > 3))
		{
		  const size_t NH_res = p1Data.residue;
		  const size_t CO_res = p2Data.residue;
		  const double inner_diameter = TPRIME::_PRIME_diameters[ 22 * p1Data.bead_type + p2Data.bead_type ];
		  const double outer_diameter = TPRIME::_PRIME_HB_aux_min_distances[3 * p1Data.bead_type + p2Data.bead_type];
		  const double bond_energy = checkTimeDependentCriteria(NH_res, CO_res, 4) ? _PRIME_HB_strength : 0;
		  return std::make_tuple(outer_diameter, inner_diameter, bond_energy, NH_res, CO_res);
		}
	      else if ((p1Data.bead_type == TPRIME::CO) && (p2Data.bead_type == TPRIME::CO) && (p1Data.location != TPRIME::CO_END) && (p2Data.location != TPRIME::CO_END))
		{
		  const size_t NH_res_1 = p2Data.residue + 1;
		  const size_t CO_res_1 = p1Data.residue;
		  const bool valid_distance_1 = (std::abs(int(NH_res_1) - int(CO_res_1)) > 3);

		  const size_t NH_res_2 = p1Data.residue + 1;
		  const size_t CO_res_2 = p2Data.residue;
		  const bool valid_distance_2 = (std::abs(int(NH_res_2) - int(CO_res_2)) > 3);
	       
		  if (valid_distance_1 || valid_distance_2)
		    {
		      const double inner_diameter = TPRIME::_PRIME_diameters[22 * p1Data.bead_type + p2Data.bead_type];
		      const double outer_diameter = TPRIME::_PRIME_HB_aux_min_distances[3 * p1Data.bead_type + p2Data.bead_type];
		      if (valid_distance_1 && checkTimeDependentCriteria(NH_res_1, CO_res_1, 3))
 			return std::make_tuple(outer_diameter, inner_diameter, _PRIME_HB_strength, NH_res_1, CO_res_1);
		      if (valid_distance_2 && checkTimeDependentCriteria(NH_res_2, CO_res_2, 3))
 			return std::make_tuple(outer_diameter, inner_diameter, _PRIME_HB_strength, NH_res_2, CO_res_2);

 		      return std::make_tuple(outer_diameter, inner_diameter, 0.0, no_HB_res, no_HB_res);
		    }
		}
	      else if ((p1Data.bead_type == TPRIME::NH) && (p2Data.bead_type == TPRIME::NH) && (p1Data.location != TPRIME::NH_END) && (p2Data.location != TPRIME::NH_END))
		{
		  //First try assuming p1 is the main NH
		  //then the main CO is pID2-1 and has a resID of p2Data.residue-1
		
		  const size_t NH_res_1 = p1Data.residue;
		  const size_t CO_res_1 = p2Data.residue - 1;
		  const bool valid_distance_1 = (std::abs(int(NH_res_1) - int(CO_res_1)) > 3);

		  const size_t NH_res_2 = p2Data.residue;
		  const size_t CO_res_2 = p1Data.residue - 1;
		  const bool valid_distance_2 = (std::abs(int(NH_res_2) - int(CO_res_2)) > 3);
	       
		  if (valid_distance_1 || valid_distance_2)
		    {
		      const double inner_diameter = TPRIME::_PRIME_diameters[22 * p1Data.bead_type + p2Data.bead_type];
		      const double outer_diameter = TPRIME::_PRIME_HB_aux_min_distances[3 * p1Data.bead_type + p2Data.bead_type];
		      if (valid_distance_1 && checkTimeDependentCriteria(NH_res_1, CO_res_1, 2))
 			return std::make_tuple(outer_diameter, inner_diameter, _PRIME_HB_strength, NH_res_1, CO_res_1);
		      if (valid_distance_2 && checkTimeDependentCriteria(NH_res_2, CO_res_2, 2))
 			return std::make_tuple(outer_diameter, inner_diameter, _PRIME_HB_strength, NH_res_2, CO_res_2);

 		      return std::make_tuple(outer_diameter, inner_diameter, 0.0, no_HB_res, no_HB_res);
		    }
		}
	      else if ((p1Data.bead_type == TPRIME::NH) && (p2Data.bead_type == TPRIME::CH) && (p1Data.location != TPRIME::NH_END) && (p2Data.location != TPRIME::CO_END) && (std::abs(int(p1Data.residue) - int(p2Data.residue)) > 3))
		{		
		  const size_t NH_res = p1Data.residue;
		  const size_t CO_res = p2Data.residue;
		  const double inner_diameter = TPRIME::_PRIME_diameters[22 * p1Data.bead_type + p2Data.bead_type];
		  const double outer_diameter = TPRIME::_PRIME_HB_aux_min_distances[3 * p1Data.bead_type + p2Data.bead_type];
		  const double bond_energy = checkTimeDependentCriteria(NH_res, CO_res, 1) ? _PRIME_HB_strength : 0;
		  return std::make_tuple(outer_diameter, inner_diameter, bond_energy, NH_res, CO_res);
		}
	      break;
	    }
          };

	//It's a standard hard-sphere interaction (no HB interactions)
	const double inner_diameter = 0.0;
	const double outer_diameter = TPRIME::_PRIME_diameters[ 22 * p1Data.bead_type + p2Data.bead_type ];
	const double bond_energy = std::numeric_limits<double>::infinity();
	return std::make_tuple(outer_diameter, inner_diameter, bond_energy, no_HB_res, no_HB_res);
      }
    else //BB-SC interaction
      {
        if (p1Data.residue == p2Data.residue) //They are [pseudo]bonded on the same residue
          {
	    const double inner_diameter = TPRIME::_PRIME_SC_BB_bond_lengths[ 22 * p1Data.bead_type + p2Data.bead_type ] * (1.0 - TPRIME::_PRIME_bond_tolerance);
	    const double outer_diameter = TPRIME::_PRIME_SC_BB_bond_lengths[ 22 * p1Data.bead_type + p2Data.bead_type ] * (1.0 + TPRIME::_PRIME_bond_tolerance);
            const double bond_energy = -std::numeric_limits<double>::infinity();
	    return std::make_tuple(outer_diameter, inner_diameter, bond_energy, no_HB_res, no_HB_res);
          }
        else
          {
            double inner_diameter = TPRIME::_PRIME_diameters[22 * p1Data.bead_type + p2Data.bead_type];
            double outer_diameter = TPRIME::_PRIME_well_diameters[22 * p1Data.bead_type + p2Data.bead_type];
            double bond_energy = TPRIME::_PRIME_well_depths[22 * p1Data.bead_type + p2Data.bead_type];

	    if (bond_energy == 0)
	      { 
		//Its a hard sphere interaction!
		bond_energy = +std::numeric_limits<double>::infinity();
		outer_diameter = inner_diameter;
		inner_diameter = 0;
	      }

            //Check for cases where it could be a "close" interaction
            if ((p2Data.residue - 1 ==  p1Data.residue) && (p1Data.bead_type == TPRIME::CO))
              { //p1 is on the residue before p2
		inner_diameter *= TPRIME::_PRIME_near_diameter_scale_factor;
		outer_diameter *= TPRIME::_PRIME_near_diameter_scale_factor;
	      }
            else if ((p2Data.residue + 1 == p1Data.residue) && (p1Data.bead_type == TPRIME::NH))
              { //p2 is on the residue after p1
		inner_diameter *= TPRIME::_PRIME_near_diameter_scale_factor;
		outer_diameter *= TPRIME::_PRIME_near_diameter_scale_factor;
              }

	    return std::make_tuple(outer_diameter, inner_diameter, bond_energy, no_HB_res, no_HB_res);
          }
      }

    M_throw() << "Missing return statement in getInteractionParameter!";
  }

  bool
  IPRIME::checkTimeDependentCriteria(const size_t NH_res, const size_t CO_res, const size_t distance_i) const
  {
    //The 5 distance criteria are labelled 0 to 4, and the current pair's applicable number is distance_i.
    //NH_ID and CO_ID give the IDs of the central NH-CO pair in the candidate hydrogen bond.

    //Check if the pair are already bonded
    if (_HBonds.count(HbondMapType::value_type(NH_res, CO_res)))
      return true; //They are!

    //Check if either pair are already within a H-Bond
    if (_HBonds.by<NH_res_ID>().count(NH_res) || _HBonds.by<CO_res_ID>().count(CO_res))
      return false; //At least one residue is already bonded, so this pair cannot form.
      
    //Search for reasons NOT to allow the bond
    if ((distance_i != 0) && !isCaptured(_topology->getBeadID(TPRIME::BeadData(TPRIME::NH, NH_res, TPRIME::MID)),
					 _topology->getBeadID(TPRIME::BeadData(TPRIME::CO, CO_res, TPRIME::MID))))
      return false;
    
    if ((distance_i != 1) && isCaptured(_topology->getBeadID(TPRIME::BeadData(TPRIME::NH, NH_res, TPRIME::MID)),
					_topology->getBeadID(TPRIME::BeadData(TPRIME::CH, CO_res, TPRIME::MID))))
      return false;

    if ((distance_i != 2) && isCaptured(_topology->getBeadID(TPRIME::BeadData(TPRIME::NH, NH_res, TPRIME::MID)),
					_topology->getBeadID(TPRIME::BeadData(TPRIME::NH, CO_res+1, TPRIME::MID))))
      return false;

    if ((distance_i != 3) && isCaptured(_topology->getBeadID(TPRIME::BeadData(TPRIME::CO, NH_res-1, TPRIME::MID)),
					_topology->getBeadID(TPRIME::BeadData(TPRIME::CO, CO_res, TPRIME::MID))))
      return false;

    if ((distance_i != 4) && isCaptured(_topology->getBeadID(TPRIME::BeadData(TPRIME::CO, NH_res, TPRIME::MID)),
					_topology->getBeadID(TPRIME::BeadData(TPRIME::CH, CO_res, TPRIME::MID))))
      return false;
    
    return true;
  }

  IntEvent
  IPRIME::getEvent(const Particle &p1, const Particle &p2) const
  {
#ifdef DYNAMO_DEBUG
    if (!Sim->dynamics->isUpToDate(p1))
      M_throw() << "Particle 1 is not up to date";

    if (!Sim->dynamics->isUpToDate(p2))
      M_throw() << "Particle 2 is not up to date";

    if (p1 == p2)
      M_throw() << "You shouldn't pass p1==p2 events to the interactions!";
#endif

    //Calculate the interaction parameters (and name them sensibly)
    const auto interaction_data = getInteractionParameters(p1.getID(), p2.getID());    
    const double outer_diameter = std::get<0>(interaction_data),
      inner_diameter = std::get<1>(interaction_data),
      bond_energy = std::get<2>(interaction_data);

    IntEvent retval(p1, p2, HUGE_VAL, NONE, *this);

    if (bond_energy == -std::numeric_limits<double>::infinity())
      { //The pair are bonded, check for events with the well edges
        double dt = Sim->dynamics->SphereSphereInRoot(p1, p2, inner_diameter);
        if (dt != HUGE_VAL) retval = IntEvent(p1, p2, dt, CORE, *this);

        dt = Sim->dynamics->SphereSphereOutRoot(p1, p2, outer_diameter);
        if (retval.getdt() > dt) retval = IntEvent(p1, p2, dt, BOUNCE, *this);
      }
    else if (bond_energy == std::numeric_limits<double>::infinity())
      { //The pair have a hard-sphere interaction, test for this event
        double dt = Sim->dynamics->SphereSphereInRoot(p1, p2, outer_diameter);
        if (dt != HUGE_VAL) retval = IntEvent(p1, p2, dt, CORE, *this);
      }
    else
      { //The pair have a square-well/shoulder interaction, test for this
	if (isCaptured(p1, p2))
	  {
	    double dt = Sim->dynamics->SphereSphereInRoot(p1, p2, inner_diameter);
	    if (dt != HUGE_VAL)
	      retval = IntEvent(p1, p2, dt, CORE, *this);
	    
	    dt = Sim->dynamics->SphereSphereOutRoot(p1, p2, outer_diameter);
	    if (retval.getdt() > dt)
	      retval = IntEvent(p1, p2, dt, STEP_OUT, *this);
	  }
	else
	  {
	    double dt = Sim->dynamics->SphereSphereInRoot(p1, p2, outer_diameter);
	    if (dt != HUGE_VAL)
	      retval = IntEvent(p1, p2, dt, STEP_IN, *this);
	  }
      }

    return retval;
  }

  PairEventData
  IPRIME::runEvent(Particle& p1, Particle& p2, const IntEvent& iEvent)
  {
    ++Sim->eventCount;

    //Calculate the interaction parameters (and name them sensibly)
    const auto interaction_data = getInteractionParameters(p1.getID(), p2.getID());    
    const double outer_diameter = std::get<0>(interaction_data),
      inner_diameter = std::get<1>(interaction_data),
      bond_energy = std::get<2>(interaction_data);
    const size_t NH_res = std::get<3>(interaction_data);
    const size_t CO_res = std::get<4>(interaction_data);
    const size_t no_HB_res = std::numeric_limits<size_t>::max();

    /*Handle the different types of interactions that can occur.*/
    PairEventData EDat;
    switch (iEvent.getType())
      {
      case CORE:
        { //CORE events occur for all types of interaction
          double coreD = inner_diameter;
	  //Its the outer_diameter for CORE events if this is a hard sphere interaction
	  if (bond_energy == std::numeric_limits<double>::infinity())
	    coreD = outer_diameter;
	  
	  EDat = Sim->dynamics->SmoothSpheresColl(iEvent, 1.0, coreD * coreD, iEvent.getType());
        }
        break;
      case BOUNCE:
        {
	  //BOUNCE events only occur for outer_diameter of a bond.
          EDat = Sim->dynamics->SmoothSpheresColl(iEvent, 1.0, outer_diameter * outer_diameter, iEvent.getType());
          break;
        }
      case STEP_IN:
	{
	  EDat = Sim->dynamics->SphereWellEvent(iEvent, -bond_energy, outer_diameter * outer_diameter, 1);
	  if (EDat.getType() != BOUNCE) {
	    //The particles have entered the well
	    if ((NH_res != no_HB_res) && (CO_res != no_HB_res) && (bond_energy != 0))
	      {
		//This is an interaction as part of the formation or breaking of a H-Bond 
		if (bond_energy < 0)
		  //Its the formation of the NH-CO bond
		  formHBond(NH_res, CO_res);
		else
		  //Its the breaking of a minimum distance criteria
		  breakHBond(NH_res, CO_res);
	      }
	    ICapture::add(p1, p2);
	  }
	  
	  break;
	}
      case STEP_OUT:
	{
	  EDat = Sim->dynamics->SphereWellEvent(iEvent, bond_energy, outer_diameter * outer_diameter, 0);
	  if (EDat.getType() != BOUNCE) {
	    //The particles are leaving the well
	    if ((NH_res != no_HB_res) && (CO_res != no_HB_res) && (bond_energy != 0))
	      {
		//This is an interaction as part of the formation or breaking of a H-Bond 
		if (bond_energy < 0)
		  //Its a NH-CO pair separating and breaking the bond
		  breakHBond(NH_res, CO_res);
		else
		  //Its a minimum distance constraint becoming satisfied, causing a H-Bond forming
		  formHBond(NH_res, CO_res);
	      }
	    ICapture::remove(p1, p2);
	  }
	  break;
	}
      default:
        M_throw() << "Unknown collision type";
      }

    return EDat;
  }

  void 
  IPRIME::formHBond(const size_t NH_res, const size_t CO_res) {
    dout << "FORMING A BOND!" << std::endl;
    auto retval = _HBonds.insert(HbondMapType::value_type(NH_res, CO_res));
    if (retval.second == false)
      M_throw() << "Failed to form a HBond";
  }

  void 
  IPRIME::breakHBond(const size_t NH_res, const size_t CO_res) {
    dout << "BREAKING A BOND!" << std::endl;
    size_t deleted_count = _HBonds.erase(HbondMapType::value_type(NH_res, CO_res));
    if (deleted_count == 0)
      M_throw() << "Failed to break a HBond";
  }

  bool 
  IPRIME::validateState(const Particle& p1, const Particle& p2, bool textoutput) const
  {
    const TPRIME::BeadData p1Data = getBeadData(p1);
    const TPRIME::BeadData p2Data = getBeadData(p2);

    //Calculate the interaction parameters (and name them sensibly)
    const auto interaction_data = getInteractionParameters(p1.getID(), p2.getID());    
    const double outer_diameter = std::get<0>(interaction_data),
      inner_diameter = std::get<1>(interaction_data),
      bond_energy = std::get<2>(interaction_data);

    if (bond_energy == -std::numeric_limits<double>::infinity())
      {//Bonded
	if (Sim->dynamics->sphereOverlap(p1, p2, inner_diameter))
	  {
	    if (textoutput)
	      derr << "Particle " << p1.getID() << " (" << TPRIME::PRIME_site_names[p1Data.bead_type] << ":"<< p1Data.residue << ") and Particle " << p2.getID()  << " (" << TPRIME::PRIME_site_names[p2Data.bead_type] << ":" << p2Data.residue << ")"
		   << " are inside the bond with an inner hard core at " << inner_diameter / Sim->units.unitLength()
		   << " but they are at a distance of " 
		   << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
		   << std::endl;
	    return true;
	  }
      
	if (!Sim->dynamics->sphereOverlap(p1, p2, outer_diameter))
	  {
	    if (textoutput)
	      derr << "Particle " << p1.getID() << " (" << TPRIME::PRIME_site_names[p1Data.bead_type] << ":"<< p1Data.residue << ") and Particle " << p2.getID()  << " (" << TPRIME::PRIME_site_names[p2Data.bead_type] << ":" << p2Data.residue << ")"
		   << " should be inside the bond with an upper limit of " << outer_diameter / Sim->units.unitLength()
		   << " but they are at a distance of " 
		   << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
		   << std::endl;
	  
	    return true;
	  }
      }
    else if (bond_energy == std::numeric_limits<double>::infinity())
      { //Hard sphere
	if (Sim->dynamics->sphereOverlap(p1, p2, outer_diameter))
	  {
	    if (textoutput)
	      derr << "Particle " << p1.getID() << " (" << TPRIME::PRIME_site_names[p1Data.bead_type] << ":"<< p1Data.residue << ") and Particle " << p2.getID()  << " (" << TPRIME::PRIME_site_names[p2Data.bead_type] << ":" << p2Data.residue << ")"
		   << " are inside the hard core at " << outer_diameter / Sim->units.unitLength()
		   << " and are at a distance of " 
		   << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
		   << std::endl;
	    return true;
	  }
      }
    else
      {
	const bool captured = isCaptured(p1, p2);

	if (captured && Sim->dynamics->sphereOverlap(p1, p2, inner_diameter))
	  {
	    if (textoutput)
	      derr << "Particle " << p1.getID() << " (" << TPRIME::PRIME_site_names[p1Data.bead_type] << ":"<< p1Data.residue << ") and Particle " << p2.getID()  << " (" << TPRIME::PRIME_site_names[p2Data.bead_type] << ":" << p2Data.residue << ")"
		   << " are inside the inner hard core of the well at " << inner_diameter / Sim->units.unitLength()
		   << " and are at a distance of " 
		   << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
		   << std::endl;
	    return true;
	  }
      
	if (captured && !Sim->dynamics->sphereOverlap(p1, p2, outer_diameter))
	  {
	    if (textoutput)
	      derr << "Particle " << p1.getID() << " (" << TPRIME::PRIME_site_names[p1Data.bead_type] << ":"<< p1Data.residue << ") and Particle " << p2.getID()  << " (" << TPRIME::PRIME_site_names[p2Data.bead_type] << ":" << p2Data.residue << ")"
		   << " are registered as being inside the well with an upper limit of " << outer_diameter / Sim->units.unitLength()
		   << " but they are at a distance of " 
		   << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
		   << std::endl;
	  
	    return true;
	  }

	if (!captured && Sim->dynamics->sphereOverlap(p1, p2, outer_diameter))
	  {
	    if (textoutput)
	      derr << "Particle " << p1.getID() << " (" << TPRIME::PRIME_site_names[p1Data.bead_type] << ":"<< p1Data.residue << ") and Particle " << p2.getID()  << " (" << TPRIME::PRIME_site_names[p2Data.bead_type] << ":" << p2Data.residue << ")"
		   << " but they are at a distance of " 
		   << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
		   << std::endl;
	  
	    return true;
	  }
      }

    return false;
  }

  void
  IPRIME::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type")  << "PRIME"
        << magnet::xml::attr("Name")  << intName
        << magnet::xml::attr("Topology") << _topology->getName()
        << magnet::xml::attr("HBStrength") << _PRIME_HB_strength
        << *range;
    
    ICapture::outputCaptureMap(XML);
    
    XML << magnet::xml::tag("HBonds");
    for (auto iter = _HBonds.left.begin(); iter != _HBonds.left.end(); ++iter)
      XML << magnet::xml::tag("Bond")
	  << magnet::xml::attr("NH")  << iter->first
	  << magnet::xml::attr("CO")  << iter->second
	  << magnet::xml::endtag("Bond");
    XML << magnet::xml::endtag("HBonds");
  }
}
