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

#include <dynamo/interactions/PRIME_BB.hpp>
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
  IPRIME_BB::IPRIME_BB(const magnet::xml::Node& XML, dynamo::Simulation* tmp):
    ICapture(tmp, NULL) //A temporary value!
  {
    operator<<(XML);
  }

  void 
  IPRIME_BB::operator<<(const magnet::xml::Node& XML)
  { 
    Interaction::operator<<(XML);

    const std::string topologyName = XML.getAttribute("Topology");
    _topology = std::dynamic_pointer_cast<TPRIME>(Sim->topology[topologyName]);
    
    if (!_topology)
      M_throw() << "For \"" << getName() << "\", Topology is not a PRIME topology.";

    ICapture::loadCaptureMap(XML);   
  }

  void 
  IPRIME_BB::initialise(size_t nID)
  {
    Interaction::initialise(nID);
    ICapture::initCaptureMap();
  }

  size_t
  IPRIME_BB::captureTest(const Particle& p1, const Particle& p2) const
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
  IPRIME_BB::getInternalEnergy() const
  { 
    //Once the capture maps are loaded just iterate through that determining energies
    double Energy = 0.0;
    for (const ICapture::value_type& IDs : *this)
      Energy += getInternalEnergy(Sim->particles[IDs.first.first], Sim->particles[IDs.first.second]);
    return Energy; 
  }

  double 
  IPRIME_BB::getInternalEnergy(const Particle& p1, const Particle& p2) const
  {
    const auto interaction_data = getInteractionParameters(p1.getID(), p2.getID());    
    const double bond_energy = std::get<2>(interaction_data);
    return bond_energy * isCaptured(p1, p2);
  }


  std::array<double, 4> 
  IPRIME_BB::getGlyphSize(size_t ID) const
  {
    return {{TPRIME::_PRIME_diameters[getBeadData(ID).first], 0, 0, 0}};
  }

  double
  IPRIME_BB::getExcludedVolume(size_t ID) const
  {
    //This calculation only includes the volumes which are always
    //excluded (i.e. the hard core)
    double diam = TPRIME::_PRIME_diameters[getBeadData(ID).first];
    return diam * diam * diam * M_PI / 6.0;
  }

  double
  IPRIME_BB::maxIntDist() const
  {
    double maxdiam = 0;
    maxdiam = std::max(maxdiam, *std::max_element(TPRIME::_PRIME_diameters, TPRIME::_PRIME_diameters + 3));
    maxdiam = std::max(maxdiam, (1.0 + TPRIME::_PRIME_bond_tolerance) * (*std::max_element(TPRIME::_PRIME_BB_bond_lengths, TPRIME::_PRIME_BB_bond_lengths + 3 * 3)));
    maxdiam = std::max(maxdiam, (1.0 + TPRIME::_PRIME_bond_tolerance) * (*std::max_element(TPRIME::_PRIME_pseudobond_lengths, TPRIME::_PRIME_pseudobond_lengths + 3 * 3)));
    maxdiam = std::max(maxdiam, (1.0 + TPRIME::_PRIME_bond_tolerance) * TPRIME::_PRIME_CH_CH_pseudobond_length);
    return maxdiam;
  }

  std::tuple<double, double, double>
  IPRIME_BB::getInteractionParameters(const size_t pID1, const size_t pID2) const
  {
    const TPRIME::BeadData p1Data = getBeadData(pID1);
    const TPRIME::BeadData p2Data = getBeadData(pID2);

    double outer_diameter; //Either the outer well diameter, or the hard-sphere diameter
    double inner_diameter; //0 if its a hard sphere
    double bond_energy = 0.0; //+inf if its a hard sphere, -inf if its a bond.

    if (p1Data.first > TPRIME::CO && p2Data.first > TPRIME::CO ) //SC-SC interaction
      {
        inner_diameter = TPRIME::_PRIME_diameters[ 22 * p1Data.first + p2Data.first ];
        outer_diameter = TPRIME::_PRIME_well_diameters[ 22 * p1Data.first + p2Data.first ];
        bond_energy    = TPRIME::_PRIME_well_depths[ 22 * p1Data.first + p2Data.first ];
      }
    else if (p1Data.first <= TPRIME::CO && p2Data.first <= TPRIME::CO) //BB-BB interaction
      {
        const size_t loc1 = p1Data.first + 3 * p1Data.second;
        const size_t loc2 = p2Data.first + 3 * p2Data.second;
        const size_t distance = std::max(loc1, loc2) - std::min(loc1, loc2);

        //This treats the special cases if they are 0,1,2, or three backbone
        //bonds apart
        switch (distance)
          {
          case 0:
            M_throw() << "Invalid backbone distance of 0";
          case 1:
            {//Every type of this interaction is a bonded interaction
              inner_diameter = TPRIME::_PRIME_BB_bond_lengths[3 * p1Data.first + p2Data.first]
                                  * (1.0 - TPRIME::_PRIME_bond_tolerance);
              outer_diameter = TPRIME::_PRIME_BB_bond_lengths[3 * p1Data.first + p2Data.first]
                                  * (1.0 + TPRIME::_PRIME_bond_tolerance);
              bond_energy    = -std::numeric_limits<double>::infinity();
            }
            break;
          case 2:
            {//Every type of this interaction is a pseudobond interaction
              inner_diameter = TPRIME::_PRIME_pseudobond_lengths[3 * p1Data.first + p2Data.first]
                                  * (1.0 - TPRIME::_PRIME_bond_tolerance);
              outer_diameter = TPRIME::_PRIME_pseudobond_lengths[3 * p1Data.first + p2Data.first]
                                  * (1.0 + TPRIME::_PRIME_bond_tolerance);
              bond_energy    = -std::numeric_limits<double>::infinity();
            }
            break;
          case 3:
            {
              //Check if this is the special pseudobond
              if ((p1Data.first == TPRIME::CH) && (p2Data.first == TPRIME::CH))
                {
                  inner_diameter = TPRIME::_PRIME_CH_CH_pseudobond_length
                                      * (1.0 - TPRIME::_PRIME_bond_tolerance);
                  outer_diameter = TPRIME::_PRIME_CH_CH_pseudobond_length
                                      * (1.0 + TPRIME::_PRIME_bond_tolerance);
                  bond_energy    = -std::numeric_limits<double>::infinity();
                }
              else
                //Close backbone-backbone hard-sphere interaction
                inner_diameter = 0.0;
                outer_diameter = TPRIME::_PRIME_diameters[ 22 * p1Data.first + p2Data.first ]
                                    * TPRIME::_PRIME_near_diameter_scale_factor;
                bond_energy = std::numeric_limits<double>::infinity();
            }
            break;
          default:
            //Backbone-backbone hard-sphere interaction
            inner_diameter = 0.0;
            outer_diameter = TPRIME::_PRIME_diameters[ 22 * p1Data.first + p2Data.first ];
            bond_energy = std::numeric_limits<double>::infinity();
            break;
          };
      }
    else //BB-SC interaction
      {
        if (p1Data.second == p2Data.second) //They are [pseudo]bonded on the same residue
          {
            if (p1Data.first <= TPRIME::CO) //p1 is BB, p2 is SC
              {
                  inner_diameter = TPRIME::_PRIME_SC_BB_bond_lengths[ 22 * p1Data.first + p2Data.first ]
                                      * (1.0 - TPRIME::_PRIME_bond_tolerance);
                  outer_diameter = TPRIME::_PRIME_SC_BB_bond_lengths[ 22 * p1Data.first + p2Data.first ]
                                      * (1.0 + TPRIME::_PRIME_bond_tolerance);
              }
            else //p2 is BB, p1 is SC
              {
                  inner_diameter = TPRIME::_PRIME_SC_BB_bond_lengths[ 22 * p2Data.first + p1Data.first ]
                                      * (1.0 - TPRIME::_PRIME_bond_tolerance);
                  outer_diameter = TPRIME::_PRIME_SC_BB_bond_lengths[ 22 * p2Data.first + p1Data.first ]
                                      * (1.0 + TPRIME::_PRIME_bond_tolerance);
              }
            bond_energy = -std::numeric_limits<double>::infinity();
          }

        else
          {
            inner_diameter = TPRIME::_PRIME_diameters[ 22 * p1Data.first + p2Data.first ];
            outer_diameter = TPRIME::_PRIME_well_diameters[ 22 * p1Data.first + p2Data.first ];
            bond_energy    = TPRIME::_PRIME_well_depths[ 22 * p1Data.first + p2Data.first ];

            //Check for cases where it could be a "close" interaction
            if (p2Data.second - 1 ==  p1Data.second) //p1 is on the residue before p2
              {
                  if ((p1Data.first > TPRIME::CO && p2Data.first == TPRIME::NH) || (p2Data.first > TPRIME::CO && p1Data.first == TPRIME::CO))
                    {
                      inner_diameter *= TPRIME::_PRIME_near_diameter_scale_factor;
                      outer_diameter *= TPRIME::_PRIME_near_diameter_scale_factor;
                    }
              }
            else if (p2Data.second + 1 == p1Data.second) //p2 is on the residue after p1
              {
                  if ((p2Data.first > TPRIME::CO && p1Data.first == TPRIME::NH) || (p1Data.first > TPRIME::CO && p2Data.first == TPRIME::CO))
                    {
                      inner_diameter *= TPRIME::_PRIME_near_diameter_scale_factor;
                      outer_diameter *= TPRIME::_PRIME_near_diameter_scale_factor;
                    }
              }
          }
      }


#ifdef DYNAMO_DEBUG
    if (bond_energy == 0.0)
      M_throw() << "Invalid bond_energy calculated, p1="<< pID1 << ", p2="<< pID2 << ", type1=" << p1Data.first << ", type2="<< p2Data.first;
#endif

    return std::make_tuple( outer_diameter, inner_diameter, bond_energy );
  }

  IntEvent
  IPRIME_BB::getEvent(const Particle &p1, const Particle &p2) const
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
  IPRIME_BB::runEvent(Particle& p1, Particle& p2, const IntEvent& iEvent)
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
  IPRIME_BB::validateState(const Particle& p1, const Particle& p2, bool textoutput) const
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
	      derr << "Particle " << p1.getID() << " (" << TPRIME::PRIME_site_names[p1Data.first] << ":"<< p1Data.second << ") and Particle " << p2.getID()  << " (" << TPRIME::PRIME_site_names[p2Data.first] << ":" << p2Data.second << ")"
		   << " are inside the bond with an inner hard core at " << inner_diameter / Sim->units.unitLength()
		   << " but they are at a distance of " 
		   << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
		   << std::endl;
	    return true;
	  }
      
	if (!Sim->dynamics->sphereOverlap(p1, p2, outer_diameter))
	  {
	    if (textoutput)
	      derr << "Particle " << p1.getID() << " (" << TPRIME::PRIME_site_names[p1Data.first] << ":"<< p1Data.second << ") and Particle " << p2.getID()  << " (" << TPRIME::PRIME_site_names[p2Data.first] << ":" << p2Data.second << ")"
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
	      derr << "Particle " << p1.getID() << " (" << TPRIME::PRIME_site_names[p1Data.first] << ":"<< p1Data.second << ") and Particle " << p2.getID()  << " (" << TPRIME::PRIME_site_names[p2Data.first] << ":" << p2Data.second << ")"
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
	      derr << "Particle " << p1.getID() << " (" << TPRIME::PRIME_site_names[p1Data.first] << ":"<< p1Data.second << ") and Particle " << p2.getID()  << " (" << TPRIME::PRIME_site_names[p2Data.first] << ":" << p2Data.second << ")"
		   << " are inside the inner hard core of the well at " << inner_diameter / Sim->units.unitLength()
		   << " and are at a distance of " 
		   << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
		   << std::endl;
	    return true;
	  }
      
	if (captured && !Sim->dynamics->sphereOverlap(p1, p2, outer_diameter))
	  {
	    if (textoutput)
	      derr << "Particle " << p1.getID() << " (" << TPRIME::PRIME_site_names[p1Data.first] << ":"<< p1Data.second << ") and Particle " << p2.getID()  << " (" << TPRIME::PRIME_site_names[p2Data.first] << ":" << p2Data.second << ")"
		   << " are registered as being inside the well with an upper limit of " << outer_diameter / Sim->units.unitLength()
		   << " but they are at a distance of " 
		   << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
		   << std::endl;
	  
	    return true;
	  }

	if (!captured && Sim->dynamics->sphereOverlap(p1, p2, outer_diameter))
	  {
	    if (textoutput)
	      derr << "Particle " << p1.getID() << " (" << TPRIME::PRIME_site_names[p1Data.first] << ":"<< p1Data.second << ") and Particle " << p2.getID()  << " (" << TPRIME::PRIME_site_names[p2Data.first] << ":" << p2Data.second << ")"
		   << " but they are at a distance of " 
		   << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
		   << std::endl;
	  
	    return true;
	  }
      }

    return false;
  }

  void
  IPRIME_BB::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type")  << "PRIME_BB"
        << magnet::xml::attr("Name")  << intName
        << magnet::xml::attr("Topology") << _topology->getName();

    ICapture::outputCaptureMap(XML);
  }
}
