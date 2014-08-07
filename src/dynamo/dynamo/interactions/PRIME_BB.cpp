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
#include <dynamo/property.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <cmath>
#include <iomanip>
#include <algorithm>

//Sources:
//[1] "α-Helix formation: Discontinuous molecular dynamics on an
//intermediate-resolution protein model", Smith and Hall (2001)
//http://onlinelibrary.wiley.com/doi/10.1002/prot.1100/full
//
//[2] "Solvent effects on the conformational transition of a model
//polyalanine peptide", Nguyen, Marchut and Hall  (2004)
//http://onlinelibrary.wiley.com/doi/10.1110/ps.04701304/full

//This is an anonymous namespace, the contents of the anonymous
//namespace will not be available outside of this file. This is used
//to place all of the settings of the model at the top of this file.
namespace {
  //const size_t NH = 0;
  //const size_t CH = 1; //Also known as C_\alpha
  //const size_t CO = 2; //Also known as C

  //All values taken from [1] unless otherwise noted.

  //These are the basic bead diameters
  const double _PRIME_diameters[] = {/*NH*/ 3.3, /*CH*/ 3.7, /*CO*/ 4.0};

  //Value taken from [2], which differs from value given in [1]
  //and produces a more realistic Ramachandran plot.
  //This is the fluctuation of the bond distance allowed.
  const double _PRIME_bond_tolerance = 0.02375;

  //This is a list of the bond lengths in the backbone. These only
  //apply to beads seperated by 1 backbone bond.
  //
  //We make this a symmetric tensor to simplify the lookups. The zero
  //entries should never be used.
  const double _PRIME_bond_lengths[] =
    /*        NH,    CH,    CO, */
    {/*NH*/0.000, 1.460, 1.330,
     /*CH*/1.460, 0.000, 1.510,
     /*CO*/1.330, 1.510, 0.000}
    ;

  //This is a list of the pseudobond lengths in the backbone. These
  //only apply to beads seperated by TWO backbone bonds.
  //
  //We make this a symmetric tensor to simplify the lookups. The zero
  //entries should never be used.
  const double _PRIME_pseudobond_lengths[] =
    /*       NH,   CH,   CO, */
    {/*NH*/0.00, 2.41, 2.45,
     /*CH*/2.41, 0.00, 2.45,
     /*CO*/2.45, 2.45, 0.00}
    ;

  //The next two constants are for the interactions between beads
  //separated by THREE backbone bonds.
  //
  //This is the special pseudobond length between the CH-CH groups. It
  //is the only pseudobond at this distance.
  const double _PRIME_CH_CH_pseudobond_length = 3.80;
  //
  //This is the scaling factor of the bead diameters if they are
  //closer than four bonds on the same chain.
  const double _PRIME_near_diameter_scale_factor = 0.75;

  //Backbone-to-backbone hydrogen bonding
  //Maximum distance for the CO and NH sites
  const double _PRIME_HB_well_diameter = 4.20;
  //
  //Minimum distance for the "auxiliary pairs"
  const double _PRIME_HB_aux_min_distances[] =
    /*       NH,   CH,   CO, */
    {/*NH*/4.74, 5.00, 0.00,
     /*CH*/5.00, 0.00, 4.86,
     /*CO*/0.00, 4.86, 4.83}
    ;
  
  enum PRIMEGroupType {
    NH, CH, CO,
    
    /*Leave this at the end */ 
    GROUP_COUNT
  };
  
  const double _PRIME_masses[] = {
    /*NH*/1.0,
    /*CH*/1.0,
    /*CO*/1.0,
  };
}

namespace dynamo {
  typedef std::unordered_map<size_t, std::pair<PRIMEGroupType, size_t> > BeadTypeMap;

  /*! \brief A class which stores the type of PRIME group that the particle corresponds to.
  
  This is a specialist class for storing the group type information.
*/
  class PRIMEGroupProperty: public Property
  {
    std::string _name;
    std::shared_ptr<BeadTypeMap> _beadTypes;
  public:
    inline PRIMEGroupProperty(const std::string name, std::shared_ptr<BeadTypeMap> map):
      Property(Units::Mass()), _name(name), _beadTypes(map) {}
    
    inline virtual const double getProperty(size_t ID) const
    { 
      
      const auto it = _beadTypes->find(ID);
      if (it == _beadTypes->end())
	M_throw() << "Do not have a PRIME bead type for particle ID " << ID;

      return _PRIME_masses[it->second.first];
    }

    inline virtual std::string getName() const 
    { return _name; }
  
    inline virtual const double getMaxValue() const
    { return *std::max_element(_PRIME_masses, _PRIME_masses + GROUP_COUNT); }

    inline virtual const double getMinValue() const 
    { return *std::min_element(_PRIME_masses, _PRIME_masses + GROUP_COUNT); }
  
    //! \sa Property::rescaleUnit
    inline virtual const void rescaleUnit(const Units::Dimension dim, const double rescale)
    {
      double factor = std::pow(rescale, _units.getUnitsPower(dim));

      if (factor && factor != 1)
	M_throw() << "Can't rescale the mass of the PRIMEGroupProperty yet!";
    }

    inline void outputParticleXMLData(magnet::xml::XmlStream& XML, const size_t pID) const
    {}
  
  
  protected:
    /*! \brief Output an XML representation of the Property to the
      passed XmlStream.
    */
    virtual void outputXML(magnet::xml::XmlStream& XML) const 
    {}
  };

  size_t
  IPRIME_BB::getType(const size_t particleID) const
  {
#ifdef DYNAMO_DEBUG
    if ((particleID < startID) || (particleID > endID))
      M_throw() << "Error, the supplied particle ID (" << particleID << " is out of range ["<< startID << "," << endID << "]";
#endif
    //The types of the ID's are constants defined above.
    //
    //This implementation will need to be improved once we generalise
    //to multiple backbones.
    return (particleID - startID) % 3;
  }

  size_t
  IPRIME_BB::getDistance(const size_t pID1, const size_t pID2) const
  {
#ifdef DYNAMO_DEBUG
    if ((pID1 < startID) || (pID2 > endID))
      M_throw() << "Error, the supplied particle ID1 (" << pID1 << " is out of range ["<< startID << "," << endID << "]";

    if ((pID2 < startID) || (pID2 > endID))
      M_throw() << "Error, the supplied particle ID2 (" << pID2 << " is out of range ["<< startID << "," << endID << "]";
#endif

    //Here, we must be careful to return only positive distances along
    //the chain. Negative values would cause the unsigned size_t to
    //roll over, plus it makes the algorithm easier to implement if
    //this is always positive.
    //
    //This implementation will need to be improved once we generalise
    //to multiple backbones.
    return std::max(pID1, pID2) - std::min(pID1, pID2);
  }


  IPRIME_BB::IPRIME_BB(const magnet::xml::Node& XML, dynamo::Simulation* tmp):
    Interaction(tmp, NULL) //A temporary value!
  {
    operator<<(XML);
  }

  void IPRIME_BB::operator<<(const magnet::xml::Node& XML)
  {
    //We don't call this operator, as we custom load our Range (it must be a linear range)
    //Interaction::operator<<(XML);
    try {
      startID = XML.getAttribute("Start").as<unsigned long>();
      endID = XML.getAttribute("End").as<unsigned long>();

      range = shared_ptr<IDPairRange>(new IDPairRangeSingle(new IDRangeRange(startID, endID)));

      intName = XML.getAttribute("Name");
    }
    catch (boost::bad_lexical_cast &) { M_throw() << "Failed a lexical cast in IPRIME_BB"; }
  }

  std::array<double, 4> 
  IPRIME_BB::getGlyphSize(size_t ID) const
  {
    return {{_PRIME_diameters[getType(ID)], 0, 0, 0}};
  }

  double
  IPRIME_BB::getExcludedVolume(size_t ID) const
  {
    //This calculation only includes the volumes which are always
    //excluded (i.e. the hard core)
    double diam = _PRIME_diameters[getType(ID)];
    return diam * diam * diam * M_PI / 6.0;
  }

  double
  IPRIME_BB::maxIntDist() const
  {
    double maxdiam = 0;
    maxdiam = std::max(maxdiam, *std::max_element(_PRIME_diameters, _PRIME_diameters + 3));
    maxdiam = std::max(maxdiam, (1.0 + _PRIME_bond_tolerance) * (*std::max_element(_PRIME_bond_lengths, _PRIME_bond_lengths + 3 * 3)));
    maxdiam = std::max(maxdiam, (1.0 + _PRIME_bond_tolerance) * (*std::max_element(_PRIME_pseudobond_lengths, _PRIME_pseudobond_lengths + 3 * 3)));
    maxdiam = std::max(maxdiam, (1.0 + _PRIME_bond_tolerance) * _PRIME_CH_CH_pseudobond_length);
    return maxdiam;
  }

  std::pair<double, bool>
  IPRIME_BB::getInteractionParameters(const size_t pID1, const size_t pID2) const
  {
    size_t p1Type = getType(pID1);
    size_t p2Type = getType(pID2);

    //We need to discover what the interaction diameter is for the
    //particles. At first, we assume its equal to the bead diameters
    double diameter = 0.5 * (_PRIME_diameters[p1Type] + _PRIME_diameters[p2Type]);
    bool bonded = false;

    //This treats the special cases if they are 0,1,2, or three backbone
    //bonds apart
    switch (getDistance(pID1, pID2))
      {
      case 0:
        M_throw() << "Invalid backbone distance of 0";
      case 1:
        {//Every type of this interaction is a bonded interaction
          diameter = _PRIME_bond_lengths[3 * p1Type + p2Type];
          bonded = true;
        }
        break;
      case 2:
        {//Every type of this interaction is a pseudobond interaction
          diameter = _PRIME_pseudobond_lengths[3 * p1Type + p2Type];
          bonded = true;
        }
        break;
      case 3:
        {
          //Check if this is the special pseudobond
          if ((p1Type == CH) && (p2Type == CH))
            {
              diameter = _PRIME_CH_CH_pseudobond_length;
              bonded = true;
            }
          else
            //Its a standard bead-bead interaction, but the diameters
            //are scaled by a factor
            diameter *= _PRIME_near_diameter_scale_factor;
        }
        break;
      default:
        //The diameter was correctly specified in the start as a
        //bead-bead interaction.
        break;
      }

#ifdef DYNAMO_DEBUG
    if (diameter == 0)
      M_throw() << "Invalid diameter calculated";
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
        double dt = Sim->dynamics->SphereSphereInRoot(p1, p2, diameter * (1.0 - _PRIME_bond_tolerance));
        if (dt != HUGE_VAL) retval = IntEvent(p1, p2, dt, CORE, *this);

        dt = Sim->dynamics->SphereSphereOutRoot(p1, p2, diameter * (1.0 + _PRIME_bond_tolerance));
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
          //a unbonded pair, the interaction distance is the diameter:
          double coreD = diameter;

          //For a bonded pair, we need to subtract the bond
          //fluctuation:
          if (bonded) coreD = diameter * (1.0 - _PRIME_bond_tolerance);

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
          double bounceD = diameter * (1.0 + _PRIME_bond_tolerance);
          EDat = Sim->dynamics->SmoothSpheresColl(iEvent, 1.0, bounceD * bounceD, iEvent.getType());
          break;
        }
      default:
        M_throw() << "Unknown collision type";
      }

    return EDat;
  }

  /*namespace{
    std::string getTypeName(size_t type)
    {
      switch (type)
	{
	case NH: return "NH";
	case CH: return "CH";
	case CO: return "CO";
	}
      M_throw() << "Invalid type";
    }
  }*/

  bool 
  IPRIME_BB::validateState(const Particle& p1, const Particle& p2, bool textoutput) const
  {
    //Calculate the interaction diameter and if the pair are bonded.
    std::pair<double, bool> interaction_data = getInteractionParameters(p1.getID(), p2.getID());
    double diameter = interaction_data.first;
    bool bonded = interaction_data.second;

    if (bonded) {
      //Inner and outer bond diameters²
      double id = diameter * (1.0 - _PRIME_bond_tolerance);
      double od = diameter * (1.0 + _PRIME_bond_tolerance);

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
        << magnet::xml::attr("Start") << startID
        << magnet::xml::attr("End")   << endID;
  }
}
