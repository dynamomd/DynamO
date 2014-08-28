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

#include <dynamo/topology/PRIME.hpp>
#include <dynamo/ranges/IDRangeRange.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <dynamo/property.hpp>
#include <dynamo/simulation.hpp>

namespace dynamo {
  /* Out-of-class definitions for the constant data arrays*/
  constexpr double TPRIME::_PRIME_diameters[];
  constexpr double TPRIME::_PRIME_masses[];  
  constexpr double TPRIME::_PRIME_BB_bond_lengths[];  
  constexpr double TPRIME::_PRIME_pseudobond_lengths[];
  constexpr double TPRIME::_PRIME_well_diameters[];
  constexpr double TPRIME::_PRIME_well_depths[];
  constexpr double TPRIME::_PRIME_SC_BB_bond_lengths[];

  const std::vector<std::string> TPRIME::PRIME_site_names {"NH", "CH", "CO", "A", "C", "D", "E", "F", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"};


  /*! \brief A class which stores the type of PRIME group that the particle corresponds to.
  
    This property is added automatically when you use the PRIME topology.
  */
  class PRIMEGroupProperty: public Property
  {
    std::string _name;
    std::shared_ptr<TPRIME::BeadTypeMap> _beadTypes;

  public:
    inline PRIMEGroupProperty(const std::string name, std::shared_ptr<TPRIME::BeadTypeMap> map):
      Property(Units::Mass()), _name(name), _beadTypes(map) {}
    
    inline virtual const double getProperty(size_t ID) const
    { 
      const auto it = _beadTypes->left.find(ID);
      if (it == _beadTypes->left.end())
	M_throw() << "Do not have a PRIME bead type for particle ID " << ID;

      return TPRIME::_PRIME_masses[it->second.first];
    }

    inline virtual std::string getName() const 
    { return _name; }
  
    inline virtual const double getMaxValue() const
    { return *std::max_element(TPRIME::_PRIME_masses, TPRIME::_PRIME_masses + TPRIME::GROUP_COUNT); }

    inline virtual const double getMinValue() const 
    { return *std::min_element(TPRIME::_PRIME_masses, TPRIME::_PRIME_masses + TPRIME::GROUP_COUNT); }
    //! \sa Property::rescaleUnit
    inline virtual const void rescaleUnit(const Units::Dimension dim, const double rescale)
    {
      const double factor = std::pow(rescale, _units.getUnitsPower(dim));

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


  TPRIME::TPRIME(const magnet::xml::Node& XML, dynamo::Simulation* Sim, unsigned int ID):
    Topology(Sim, ID),
    _types(new BeadTypeMap)
  {
    TPRIME::operator<<(XML);
    Sim->_properties.addNamedProperty(std::shared_ptr<Property>(new PRIMEGroupProperty(_name, _types)));
  }

  void 
  TPRIME::operator<<(const magnet::xml::Node& XML) 
  {
    Topology::operator<<(XML);

    size_t residue = 0;
    for (magnet::xml::Node node = XML.findNode("Molecule"); node.valid(); ++node)
      {
	//Store the loaded data for output later (the internal
	//representation is not easy to convert back)
	_configData.push_back(std::make_pair(node.getAttribute("StartID").as<size_t>(),
					     node.getAttribute("Sequence").as<std::string>()));

	//Create an internal representation which allows fast look-ups
	const size_t startID = node.getAttribute("StartID").as<size_t>();
	size_t ID = startID;
	for (char letter: node.getAttribute("Sequence").as<std::string>())
	  {
	    _types->insert(BeadTypeMap::value_type(ID++, BeadData(NH, residue)));
	    _types->insert(BeadTypeMap::value_type(ID++, BeadData(CH, residue)));
	    _types->insert(BeadTypeMap::value_type(ID++, BeadData(CO, residue)));
	    
	    switch (letter) {
	    case 'A': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::A, residue))); break;
	    case 'C': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::C, residue))); break;
	    case 'D': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::D, residue))); break;
	    case 'E': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::E, residue))); break;
	    case 'F': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::F, residue))); break;
	    case 'G': /*This residue has no side chain*/ break;
	    case 'H': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::H, residue))); break;
	    case 'I': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::I, residue))); break;
	    case 'K': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::K, residue))); break;
	    case 'L': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::L, residue))); break;
	    case 'M': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::M, residue))); break;
	    case 'N': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::N, residue))); break;
	    case 'P': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::P, residue))); break;
	    case 'Q': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::Q, residue))); break;
	    case 'R': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::R, residue))); break;
	    case 'S': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::S, residue))); break;
	    case 'T': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::T, residue))); break;
	    case 'V': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::V, residue))); break;
	    case 'W': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::W, residue))); break;
	    case 'Y': _types->insert(BeadTypeMap::value_type(ID++, BeadData(TPRIME::Y, residue))); break;
	    default:
	      M_throw() << "Unrecognised PRIME group type " << letter;
	    }
	    ++residue;
	  }

	ranges.push_back(shared_ptr<IDRange>(new IDRangeRange(startID, ID)));
      }
  }

  void 
  TPRIME::outputXML(magnet::xml::XmlStream& XML) const 
  {
    XML << magnet::xml::attr("Name") << _name;
    XML << magnet::xml::attr("Type") << "PRIME";

    for (const auto& entry: _configData)
      XML << magnet::xml::tag("Molecule")
	  << magnet::xml::attr("StartID") << entry.first
	  << magnet::xml::attr("Sequence") << entry.second
	  << magnet::xml::endtag("Molecule");
  }
}
