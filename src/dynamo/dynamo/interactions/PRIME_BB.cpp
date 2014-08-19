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
    Interaction(tmp, NULL) //A temporary value!
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
  IPRIME_BB::getInteractionParameters2(const size_t pID1, const size_t pID2) const
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
        //Backbone-backbone interaction!
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
                outer_diameter = 0.5 * (TPRIME::_PRIME_diameters[p1Data.first] + TPRIME::_PRIME_diameters[p2Data.first]);
                outer_diameter *= TPRIME::_PRIME_near_diameter_scale_factor;
                bond_energy = std::numeric_limits<double>::infinity();
            }
            break;
          default:
            //Backbone-backbone hard-sphere interaction
            inner_diameter = 0.0;
            outer_diameter = 0.5 * (TPRIME::_PRIME_diameters[p1Data.first] + TPRIME::_PRIME_diameters[p2Data.first]);
            bond_energy = std::numeric_limits<double>::infinity();
            break;
          }
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
            if (p2Data.second - p1Data.second == 1) //p1 is on the residuce before p2
              {
                  if ((p1Data.first > TPRIME::CO && p2Data.first == TPRIME::NH) || (p2Data.first > TPRIME::CO && p1Data.first == TPRIME::CO))
                    {
                      inner_diameter *= TPRIME::_PRIME_near_diameter_scale_factor;
                      outer_diameter *= TPRIME::_PRIME_near_diameter_scale_factor;
                    }
              }
            else if (p2Data.second - p1Data.second == -1) //p2 is on the residue before p1
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
              M_throw() << "Invalid diameter calculated, p1="<< pID1 << ", p2="<<pID2 << ", distance="<<distance << ", type1=" << p1Data.first << ", type2="<< p2Data.first;
#endif

    return std::make_tuple( outer_diameter, inner_diameter, bond_energy );
  }

  std::pair<double, bool>
  IPRIME_BB::getInteractionParameters(const size_t pID1, const size_t pID2) const
  {
    const TPRIME::BeadData p1Type = getBeadData(pID1);
    const TPRIME::BeadData p2Type = getBeadData(pID2);

#ifdef DYNAMO_DEBUG
    if ((p1Type.first > TPRIME::CO) || (p2Type.first > TPRIME::CO))
      M_throw() << "Error! Can't do backbone-sidechain or sidechain-sidechain this way yet!";
#endif

    //We need to discover what the interaction diameter is for the
    //particles. At first, we assume its equal to the bead diameters
    double diameter = 0.5 * (TPRIME::_PRIME_diameters[p1Type.first] + TPRIME::_PRIME_diameters[p2Type.first]);
    bool bonded = false;

    //For backbone atoms only! We can calculate the distance like so:
    const size_t loc1 = p1Type.first + 3 * p1Type.second;
    const size_t loc2 = p2Type.first + 3 * p2Type.second;
    const size_t distance = std::max(loc1, loc2) - std::min(loc1, loc2);

    //This treats the special cases if they are 0,1,2, or three backbone
    //bonds apart
    switch (distance)
      {
      case 0:
        M_throw() << "Invalid backbone distance of 0";
      case 1:
        {//Every type of this interaction is a bonded interaction
          diameter = TPRIME::_PRIME_BB_bond_lengths[3 * p1Type.first + p2Type.first];
          bonded = true;
        }
        break;
      case 2:
        {//Every type of this interaction is a pseudobond interaction
          diameter = TPRIME::_PRIME_pseudobond_lengths[3 * p1Type.first + p2Type.first];
          bonded = true;
        }
        break;
      case 3:
        {
          //Check if this is the special pseudobond
          if ((p1Type.first == TPRIME::CH) && (p2Type.first == TPRIME::CH))
            {
              diameter = TPRIME::_PRIME_CH_CH_pseudobond_length;
              bonded = true;
            }
          else
            //Its a standard bead-bead interaction, but the diameters
            //are scaled by a factor
            diameter *= TPRIME::_PRIME_near_diameter_scale_factor;
        }
        break;
      default:
        //The diameter was correctly specified in the start as a
        //bead-bead interaction.
        break;
      }

#ifdef DYNAMO_DEBUG
    if (diameter == 0)
      M_throw() << "Invalid diameter calculated, p1="<< pID1 << ", p2="<<pID2 << ", distance="<<distance << ", type1=" << p1Type.first << ", type2="<< p2Type.first;
#endif

    return std::make_pair(diameter, bonded);
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

    //Calculate the interaction diameter and if the pair are bonded
    std::pair<double, bool> interaction_data = getInteractionParameters(p1.getID(), p2.getID());
    double diameter = interaction_data.first;
    bool bonded = interaction_data.second;

    //Now that the diameters are known, and we know if they are bonded
    //or not, we can do event detection. If it is bonded
    //(bonded=true), the inner and outer interactions are at:
    //diameter*(1+-_PRIME_bond_tolerance). If not, its a hard core
    //interaction at the diameter.

    IntEvent retval(p1, p2, HUGE_VAL, NONE, *this);
    if (bonded)
      {
        double dt = Sim->dynamics->SphereSphereInRoot(p1, p2, diameter * (1.0 - TPRIME::_PRIME_bond_tolerance));
        if (dt != HUGE_VAL) retval = IntEvent(p1, p2, dt, CORE, *this);

        dt = Sim->dynamics->SphereSphereOutRoot(p1, p2, diameter * (1.0 + TPRIME::_PRIME_bond_tolerance));
        if (retval.getdt() > dt) retval = IntEvent(p1, p2, dt, BOUNCE, *this);
      }
    else
      {
        double dt = Sim->dynamics->SphereSphereInRoot(p1, p2, diameter);
        if (dt != HUGE_VAL) retval = IntEvent(p1, p2, dt, CORE, *this);
      }

    return retval;
  }

  PairEventData
  IPRIME_BB::runEvent(Particle& p1, Particle& p2, const IntEvent& iEvent)
  {
    ++Sim->eventCount;

    //Calculate the interaction diameter and if the pair are bonded.
    const std::pair<double, bool> interaction_data = getInteractionParameters(p1.getID(), p2.getID());
    const double diameter = interaction_data.first;
    const bool bonded = interaction_data.second;

    /*Handle the different types of interactions that can occur.*/
    PairEventData EDat;
    switch (iEvent.getType())
      {
      case CORE:
        { //This is either a core for a unbonded or a bonded pair. For
          //an unbonded pair, the interaction distance is the diameter:
          double coreD = diameter;

          //For a bonded pair, we need to subtract the bond
          //fluctuation:
          if (bonded) coreD = diameter * (1.0 - TPRIME::_PRIME_bond_tolerance);

          //Finally, run the event
          EDat = Sim->dynamics->SmoothSpheresColl(iEvent, 1.0, coreD * coreD, iEvent.getType());
        }
        break;
      case BOUNCE:
        {
#ifdef DYNAMO_DEBUG
          if (!bonded) M_throw() << "Only bonded particles can undergo a BOUNCE event";
#endif
          //As this is the outer bond event, we know its diameter to
          //be:
          double bounceD = diameter * (1.0 + TPRIME::_PRIME_bond_tolerance);
          EDat = Sim->dynamics->SmoothSpheresColl(iEvent, 1.0, bounceD * bounceD, iEvent.getType());
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
    //Calculate the interaction diameter and if the pair are bonded.
    std::pair<double, bool> interaction_data = getInteractionParameters(p1.getID(), p2.getID());
    double diameter = interaction_data.first;
    bool bonded = interaction_data.second;

    if (bonded) {
      //Inner and outer bond diameters²
      double id = diameter * (1.0 - TPRIME::_PRIME_bond_tolerance);
      double od = diameter * (1.0 + TPRIME::_PRIME_bond_tolerance);

      if (Sim->dynamics->sphereOverlap(p1, p2, id))
	{
	  if (textoutput)
	    derr << "Particle " << p1.getID() << " and Particle " << p2.getID() 
		 << " are inside the bond with an inner hard core at " << id / Sim->units.unitLength()
		 << " but they are at a distance of " 
		 << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
		 << std::endl;
	  return true;
	}
      
      if (!Sim->dynamics->sphereOverlap(p1, p2, od))
	{
	  if (textoutput)
	    derr << "Particle " << p1.getID() << " and Particle " << p2.getID() 
		 << " registered as being inside the bond with an upper limit of " << od / Sim->units.unitLength()
		 << " but they are at a distance of " 
		 << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
		 << std::endl;
	  
	  return true;
	}
    }
    else
      if (Sim->dynamics->sphereOverlap(p1, p2, diameter))
	{
	  if (textoutput)
	    derr << "Particle " << p1.getID() << " and Particle " << p2.getID() 
		 << " are inside the hard core at " << diameter / Sim->units.unitLength()
		 << " and are at a distance of " 
		 << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
		 << std::endl;
	  return true;
	}

    return false;
  }

  void
  IPRIME_BB::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type")  << "PRIME_BB"
        << magnet::xml::attr("Name")  << intName
        << magnet::xml::attr("Topology") << _topology->getName();
  }
}
