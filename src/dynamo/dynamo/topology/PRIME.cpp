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
      const auto it = _beadTypes->find(ID);
      if (it == _beadTypes->end())
	M_throw() << "Do not have a PRIME bead type for particle ID " << ID;

      return TPRIME::_PRIME_masses[it->second.second];
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

    for (magnet::xml::Node node = XML.findNode("Molecule"); node.valid(); ++node)
      {
	const std::string seq = node.getAttribute("Sequence").as<std::string>();	
	const size_t startID = node.getAttribute("StartID").as<size_t>();
	size_t endID = startID;
	size_t residue = 0;
	for (char letter: seq)
	  {
	    (*_types)[endID+0] = BeadData(NH, residue);
	    (*_types)[endID+1] = BeadData(CH, residue);
	    (*_types)[endID+2] = BeadData(CO, residue);
	    
	    if (letter == 'G')
	      endID += 3;
	    else
	      {
		switch (letter) {
		case 'A': (*_types)[endID+3] = BeadData(A, residue); break;
		case 'C': (*_types)[endID+3] = BeadData(C, residue); break;
		case 'D': (*_types)[endID+3] = BeadData(D, residue); break;
		case 'E': (*_types)[endID+3] = BeadData(E, residue); break;
		case 'F': (*_types)[endID+3] = BeadData(F, residue); break;
		case 'H': (*_types)[endID+3] = BeadData(H, residue); break;
		case 'I': (*_types)[endID+3] = BeadData(I, residue); break;
		case 'K': (*_types)[endID+3] = BeadData(K, residue); break;
		case 'L': (*_types)[endID+3] = BeadData(L, residue); break;
		case 'M': (*_types)[endID+3] = BeadData(M, residue); break;
		case 'N': (*_types)[endID+3] = BeadData(N, residue); break;
		case 'P': (*_types)[endID+3] = BeadData(P, residue); break;
		case 'Q': (*_types)[endID+3] = BeadData(Q, residue); break;
		case 'R': (*_types)[endID+3] = BeadData(R, residue); break;
		case 'S': (*_types)[endID+3] = BeadData(S, residue); break;
		case 'T': (*_types)[endID+3] = BeadData(T, residue); break;
		case 'V': (*_types)[endID+3] = BeadData(V, residue); break;
		case 'W': (*_types)[endID+3] = BeadData(W, residue); break;
		case 'Y': (*_types)[endID+3] = BeadData(Y, residue); break;
		default:
		  M_throw() << "Unrecognised PRIME group type " << letter;
		}
		endID += 4;
	      }
	    ++residue;
	  }

	ranges.push_back(shared_ptr<IDRange>(new IDRangeRange(startID, endID)));
      }
  }

  void 
  TPRIME::outputXML(magnet::xml::XmlStream& XML) const 
  {
    XML << magnet::xml::attr("Name") << _name;
    XML << magnet::xml::attr("Type") << "PRIME";

    for (const shared_ptr<IDRange>& plugPtr : ranges)
      XML << magnet::xml::tag("Molecule") << plugPtr
	  << magnet::xml::endtag("Molecule");
  }
}
