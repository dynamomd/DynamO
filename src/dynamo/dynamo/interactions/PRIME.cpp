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
  }

  void 
  IPRIME::initialise(size_t nID)
  {
    Interaction::initialise(nID);
    ICapture::initCaptureMap();
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

    double outer_diameter; //Either the outer well diameter, or the hard-sphere diameter
    double inner_diameter; //0 if its a hard sphere
    double bond_energy = 0.0; //+inf if its a hard sphere, -inf if its a bond.

    //ResIDs of the main NH-CO pair for HB interactions
    size_t NH_res = 0;
    size_t CO_res = 0;

    if (p1Data.bead_type > TPRIME::CO) //SC-SC interaction (as p2Data.bead_type >= p1Data.bead_type)
      {
        inner_diameter = TPRIME::_PRIME_diameters[ 22 * p1Data.bead_type + p2Data.bead_type];
        outer_diameter = TPRIME::_PRIME_well_diameters[ 22 * p1Data.bead_type + p2Data.bead_type];
        bond_energy    = TPRIME::_PRIME_well_depths[ 22 * p1Data.bead_type + p2Data.bead_type];
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
            {//Every type of this interaction is a bonded interaction
              inner_diameter = TPRIME::_PRIME_BB_bond_lengths[3 * p1Data.bead_type + p2Data.bead_type]
                                  * (1.0 - TPRIME::_PRIME_bond_tolerance);
              outer_diameter = TPRIME::_PRIME_BB_bond_lengths[3 * p1Data.bead_type + p2Data.bead_type]
                                  * (1.0 + TPRIME::_PRIME_bond_tolerance);
              bond_energy    = -std::numeric_limits<double>::infinity();
            }
            break;
          case 2:
            {//Every type of this interaction is a pseudobond interaction
              inner_diameter = TPRIME::_PRIME_pseudobond_lengths[3 * p1Data.bead_type + p2Data.bead_type]
                                  * (1.0 - TPRIME::_PRIME_bond_tolerance);
              outer_diameter = TPRIME::_PRIME_pseudobond_lengths[3 * p1Data.bead_type + p2Data.bead_type]
                                  * (1.0 + TPRIME::_PRIME_bond_tolerance);
              bond_energy    = -std::numeric_limits<double>::infinity();
            }
            break;
          case 3:
            {
              //Check if this is the special pseudobond
              if ((p1Data.bead_type == TPRIME::CH) && (p2Data.bead_type == TPRIME::CH))
                {
                  inner_diameter = TPRIME::_PRIME_CH_CH_pseudobond_length
                                      * (1.0 - TPRIME::_PRIME_bond_tolerance);
                  outer_diameter = TPRIME::_PRIME_CH_CH_pseudobond_length
                                      * (1.0 + TPRIME::_PRIME_bond_tolerance);
                  bond_energy    = -std::numeric_limits<double>::infinity();
                }
              else
		{
		  //Close backbone-backbone hard-sphere interaction
		  inner_diameter = 0.0;
		  outer_diameter = TPRIME::_PRIME_diameters[ 22 * p1Data.bead_type + p2Data.bead_type ]
		    * TPRIME::_PRIME_near_diameter_scale_factor;
		  bond_energy = std::numeric_limits<double>::infinity();
		}
	    }
            break;
          default:
            //Backbone-backbone hard-sphere or hydrogen bond interaction

            //Split the criteria based on time-dependence. If the time-independent criteria are met,
            //we want to track the pair's capture state and evaluate the time-dependent criteria.
            //If the time-dependent criteria are met, we want to turn the energy of the interaction on.
            bool timeIndependentHBCriteria = false;
            bool timeDependentHBCriteria = false;

            inner_diameter = TPRIME::_PRIME_diameters[ 22 * p1Data.bead_type + p2Data.bead_type ];

            //Determine if HB
            //If the time-independent criteria are met, it's not a hard-sphere and we track pairs.
            //If the time-dependent criteria are met, the bond_energy is nonzero
            if (p1Data.bead_type == TPRIME::CO)
              {
                if (p2Data.bead_type == TPRIME::CO)
                  {
                    outer_diameter = TPRIME::_PRIME_HB_aux_min_distances[3 * p1Data.bead_type + p2Data.bead_type];

                    //First try assuming pID1 is the main CO
                    //then the main NH is pID2+1 and has a resID of p2Data.residue+1
                    NH_res = p2Data.residue + 1;
                    CO_res = p1Data.residue;

                    if (p1Data.location != TPRIME::CO_END && p2Data.location != TPRIME::CO_END)
                      {
                        if (abs( NH_res - CO_res ) > 3)
                          {
                            timeIndependentHBCriteria = true;
                            timeDependentHBCriteria = checkTimeDependentCriteria(NH_res, CO_res, 3);
                          }

                        //If time dependent criteria not met, try assuming pID2 is the main CO
                        //then the main NH is pID1+1 and has a resID of p1Data.residue+1
                        if (!timeDependentHBCriteria)
                          {
                            NH_res = p1Data.residue + 1;
                            CO_res = p2Data.residue;

                            if (abs( NH_res - CO_res ) > 3)
                              {
                                timeIndependentHBCriteria = true;

                                timeDependentHBCriteria = checkTimeDependentCriteria(NH_res, CO_res, 3);
                              }
                          }

                        if (timeDependentHBCriteria)
                          {
                            bond_energy = _PRIME_HB_strength;
                          }
                      }
                  }
                else if (p2Data.bead_type == TPRIME::NH)
                  {
                    // this is the primary NH-CO pair.
                    outer_diameter = TPRIME::_PRIME_HB_well_diameter;

                    NH_res = p2Data.residue;
                    CO_res = p1Data.residue;

                    if ( p1Data.location != TPRIME::CO_END && p2Data.location != TPRIME::NH_END &&
                         (abs( NH_res - CO_res ) > 3) )
                      {
                        timeIndependentHBCriteria = true;
                        timeDependentHBCriteria = checkTimeDependentCriteria(NH_res, CO_res, 0);
                      }

                    if (timeDependentHBCriteria)
                      {
                        bond_energy = -_PRIME_HB_strength;
                      }

                  }
                else if (p2Data.bead_type == TPRIME::CH)
                  {
                    outer_diameter = TPRIME::_PRIME_HB_aux_min_distances[3 * p1Data.bead_type + p2Data.bead_type];

                    NH_res = p2Data.residue;
                    CO_res = p1Data.residue;

                    if ( p1Data.location != TPRIME::CO_END && p2Data.location != TPRIME::NH_END &&
                        (abs( NH_res - CO_res ) > 3) )
                      {
                        timeIndependentHBCriteria = true;
                        timeDependentHBCriteria = checkTimeDependentCriteria(NH_res, CO_res, 4);
                      }

                    if (timeDependentHBCriteria)
                      {
                        bond_energy = _PRIME_HB_strength;
                      }
                  }
              }
            else if (p2Data.bead_type == TPRIME::CO)
              {
                if (p1Data.bead_type == TPRIME::NH)
                  {
                    // this is the primary NH-CO pair.
                    outer_diameter = TPRIME::_PRIME_HB_well_diameter;

                    NH_res = p1Data.residue;
                    CO_res = p2Data.residue;

                    if ( p1Data.location != TPRIME::CO_END && p2Data.location != TPRIME::NH_END &&
                         (abs( NH_res - CO_res ) > 3) )
                      {
                        timeIndependentHBCriteria = true;
                        timeDependentHBCriteria = checkTimeDependentCriteria(NH_res, CO_res, 0);
                      }

                    if (timeDependentHBCriteria)
                      {
                        bond_energy = -_PRIME_HB_strength;
                      }

                  }
                else if (p1Data.bead_type == TPRIME::CH)
                  {
                    outer_diameter = TPRIME::_PRIME_HB_aux_min_distances[3 * p1Data.bead_type + p2Data.bead_type];

                    NH_res = p1Data.residue;
                    CO_res = p2Data.residue;

                    if ( p2Data.location != TPRIME::CO_END && p1Data.location != TPRIME::NH_END &&
                        (abs( NH_res - CO_res ) > 3) )
                      {
                        timeIndependentHBCriteria = true;
                        timeDependentHBCriteria = checkTimeDependentCriteria(NH_res, CO_res, 4);
                      }

                    if (timeDependentHBCriteria)
                      {
                        bond_energy = _PRIME_HB_strength;
                      }
                  }
              }
            else if (p1Data.bead_type == TPRIME::NH)
              {
                if (p2Data.bead_type == TPRIME::NH)
                  {
                    outer_diameter = TPRIME::_PRIME_HB_aux_min_distances[3 * p1Data.bead_type + p2Data.bead_type];

                    //First try assuming p1 is the main NH
                    //then the main CO is pID2-1 and has a resID of p2Data.residue-1
                    NH_res = p1Data.residue;
                    CO_res = p2Data.residue - 1;

                    if (p1Data.location != TPRIME::NH_END && p2Data.location != TPRIME::NH_END)
                      {
                        if (abs( NH_res - CO_res ) > 3)
                          {
                            timeIndependentHBCriteria = true;
                            timeDependentHBCriteria = checkTimeDependentCriteria(NH_res, CO_res, 2);
                          }

                        //If time dependent criteria not met, try assuming p2 is the main NH
                        //then the main CO is pID1-1 and has a resID of p1Data.residue-1
                        if (!timeDependentHBCriteria)
                          {

                            NH_res = p2Data.residue;
                            CO_res = p1Data.residue - 1;

                            if (abs( NH_res - CO_res ) > 3)
                              {
                                timeIndependentHBCriteria = true;
                                timeDependentHBCriteria = checkTimeDependentCriteria(NH_res, CO_res, 2);
                              }
                          }

                        if (timeDependentHBCriteria)
                          {
                            bond_energy = _PRIME_HB_strength;
                          }
                      }
                  }
                else if (p2Data.bead_type == TPRIME::CH)
                  {
                    outer_diameter = TPRIME::_PRIME_HB_aux_min_distances[3 * p1Data.bead_type + p2Data.bead_type];

                    NH_res = p1Data.residue;
                    CO_res = p2Data.residue;

                    //pID1 is the main CO. The main NH is pID2-1
                    if ( p1Data.location != TPRIME::CO_END && p2Data.location != TPRIME::NH_END &&
                        (abs( NH_res - CO_res ) > 3) )
                      {
                        timeIndependentHBCriteria = true;
                        timeDependentHBCriteria = checkTimeDependentCriteria(NH_res, CO_res, 1);
                      }

                    if (timeDependentHBCriteria)
                      {
                        bond_energy = _PRIME_HB_strength;
                      }
                  }
              }
            else if (p2Data.bead_type == TPRIME::NH)
              {
                //p1 must be CH, all other possibilities exhausted
                outer_diameter = TPRIME::_PRIME_HB_aux_min_distances[3 * p1Data.bead_type + p2Data.bead_type];

                NH_res = p2Data.residue;
                CO_res = p1Data.residue;

                if ( p1Data.location != TPRIME::CO_END && p2Data.location != TPRIME::NH_END &&
                    (abs( NH_res - CO_res ) > 3) )
                  {
                    timeIndependentHBCriteria = true;
                    timeDependentHBCriteria = checkTimeDependentCriteria(NH_res, CO_res, 1);
                  }

                if (timeDependentHBCriteria)
                  {
                    bond_energy = _PRIME_HB_strength;
                  }
              }

            if (!timeIndependentHBCriteria)
              {
                //It's a hard-sphere interaction
                inner_diameter = 0.0;
                outer_diameter = TPRIME::_PRIME_diameters[ 22 * p1Data.bead_type + p2Data.bead_type ];
                bond_energy = std::numeric_limits<double>::infinity();
              }

            break;
          };
      }
    else //BB-SC interaction
      {
        if (p1Data.residue == p2Data.residue) //They are [pseudo]bonded on the same residue
          {
            if (p1Data.bead_type <= TPRIME::CO) //p1 is BB, p2 is SC
              {
                  inner_diameter = TPRIME::_PRIME_SC_BB_bond_lengths[ 22 * p1Data.bead_type + p2Data.bead_type ]
                                      * (1.0 - TPRIME::_PRIME_bond_tolerance);
                  outer_diameter = TPRIME::_PRIME_SC_BB_bond_lengths[ 22 * p1Data.bead_type + p2Data.bead_type ]
                                      * (1.0 + TPRIME::_PRIME_bond_tolerance);
              }
            else //p2 is BB, p1 is SC
              {
                  inner_diameter = TPRIME::_PRIME_SC_BB_bond_lengths[ 22 * p2Data.bead_type + p1Data.bead_type ]
                                      * (1.0 - TPRIME::_PRIME_bond_tolerance);
                  outer_diameter = TPRIME::_PRIME_SC_BB_bond_lengths[ 22 * p2Data.bead_type + p1Data.bead_type ]
                                      * (1.0 + TPRIME::_PRIME_bond_tolerance);
              }
            bond_energy = -std::numeric_limits<double>::infinity();
          }
        else
          {
            inner_diameter = TPRIME::_PRIME_diameters[22 * p1Data.bead_type + p2Data.bead_type];
            outer_diameter = TPRIME::_PRIME_well_diameters[22 * p1Data.bead_type + p2Data.bead_type];
            bond_energy    = TPRIME::_PRIME_well_depths[22 * p1Data.bead_type + p2Data.bead_type];

	    if (bond_energy == 0)
	      { 
		//Its a hard sphere interaction!
		bond_energy = +std::numeric_limits<double>::infinity();
		outer_diameter = inner_diameter;
		inner_diameter = 0;
	      }

            //Check for cases where it could be a "close" interaction
            if (p2Data.residue - 1 ==  p1Data.residue) //p1 is on the residue before p2
              {
                  if ((p1Data.bead_type > TPRIME::CO && p2Data.bead_type == TPRIME::NH) || (p2Data.bead_type > TPRIME::CO && p1Data.bead_type == TPRIME::CO))
                    {
                      inner_diameter *= TPRIME::_PRIME_near_diameter_scale_factor;
                      outer_diameter *= TPRIME::_PRIME_near_diameter_scale_factor;
                    }
              }
            else if (p2Data.residue + 1 == p1Data.residue) //p2 is on the residue after p1
              {
                  if ((p2Data.bead_type > TPRIME::CO && p1Data.bead_type == TPRIME::NH) || (p1Data.bead_type > TPRIME::CO && p2Data.bead_type == TPRIME::CO))
                    {
                      inner_diameter *= TPRIME::_PRIME_near_diameter_scale_factor;
                      outer_diameter *= TPRIME::_PRIME_near_diameter_scale_factor;
                    }
              }
          }
      }


#ifdef DYNAMO_DEBUG
    if (bond_energy == 0.0)
      M_throw() << "Invalid bond_energy calculated, p1="<< pID1 << ", p2="<< pID2 << ", type1=" << p1Data.bead_type << ", type2="<< p2Data.bead_type;
#endif

    return std::make_tuple( outer_diameter, inner_diameter, bond_energy, NH_res, CO_res);
  }

  bool
  IPRIME::checkTimeDependentCriteria(const size_t NH_res, const size_t CO_res, const size_t distance_i) const
  {
    //The 5 distance criteria are labelled 0 to 4, and the current pair's applicable number is distance_i.
    //NH_ID and CO_ID give the IDs of the central NH-CO pair in the candidate hydrogen bond.

    bool satisfied = ( _HBondedNHs.find(NH_res) == _HBondedNHs.end() &&
                       _HBondedCOs.find(CO_res) == _HBondedCOs.end() );

    //NH-CO
    if (satisfied && distance_i != 0)
      {
        satisfied = isCaptured( _topology->getBeadID(TPRIME::BeadData(TPRIME::NH, NH_res, TPRIME::MID)),
                                _topology->getBeadID(TPRIME::BeadData(TPRIME::CO, CO_res, TPRIME::MID)) );
      }
    //NH-CH
    if (satisfied && distance_i != 1)
      {
        satisfied = !isCaptured( _topology->getBeadID(TPRIME::BeadData(TPRIME::NH, NH_res, TPRIME::MID)),
                                _topology->getBeadID(TPRIME::BeadData(TPRIME::CH, CO_res, TPRIME::MID)) );
      }
    //NH-NH
    if (satisfied && distance_i != 2)
      {
        satisfied = !isCaptured( _topology->getBeadID(TPRIME::BeadData(TPRIME::NH, NH_res, TPRIME::MID)),
                                _topology->getBeadID(TPRIME::BeadData(TPRIME::NH, CO_res+1, TPRIME::MID)) );
      }
    //CO-CO
    if (satisfied && distance_i != 3)
      {
        satisfied = !isCaptured( _topology->getBeadID(TPRIME::BeadData(TPRIME::CO, NH_res-1, TPRIME::MID)),
                                _topology->getBeadID(TPRIME::BeadData(TPRIME::CO, CO_res, TPRIME::MID)) );
      }
    //CO-CH
    if (satisfied && distance_i != 4)
      {
        satisfied = !isCaptured( _topology->getBeadID(TPRIME::BeadData(TPRIME::CH, NH_res, TPRIME::MID)),
                                _topology->getBeadID(TPRIME::BeadData(TPRIME::CO, CO_res, TPRIME::MID)) );
      }

    return satisfied;
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
	  if (EDat.getType() != BOUNCE) ICapture::add(p1, p2);
	  break;
	}
      case STEP_OUT:
	{
	  EDat = Sim->dynamics->SphereWellEvent(iEvent, bond_energy, outer_diameter * outer_diameter, 0);
	  if (EDat.getType() != BOUNCE) ICapture::remove(p1, p2);
	  break;
	}
      default:
        M_throw() << "Unknown collision type";
      }

    return EDat;
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
  }
}
