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

#pragma once

#include <dynamo/base.hpp>
#include <dynamo/ranges/IDRange.hpp>
#include <dynamo/simulation.hpp>
#include <string>

namespace magnet { namespace xml { class Node; } }
namespace xml { class XmlStream; }

namespace dynamo {
  class Particle;
  class Interaction;
  class RenderObj;

  class Species: public dynamo::SimBase
  {
  public:
    virtual ~Species();

    inline bool isSpecies(const Particle& p1) const { return range->isInRange(p1); }  
    inline const double getMass(size_t ID) const { return _mass->getProperty(ID); }
    inline unsigned long getCount() const { return range->size(); }
    inline unsigned int getID() const { return ID; }
    inline const std::string& getName() const { return spName; }
    inline const shared_ptr<IDRange>& getRange() const { return range; }
    virtual double getScalarMomentOfInertia(size_t ID) const = 0;

    virtual void operator<<(const magnet::xml::Node&) = 0;

    virtual void initialise() = 0;

    friend magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream&, const Species&);
  
    static shared_ptr<Species> getClass(const magnet::xml::Node&, dynamo::Simulation*, size_t);

  protected:
    template<class T1>
    Species(dynamo::Simulation* tmp, std::string name, IDRange* nr, T1 mass, std::string nName, unsigned int nID):
      SimBase(tmp, name),
      _mass(Sim->_properties.getProperty(mass, Property::Units::Mass())),
      range(nr),
      spName(nName),
      ID(nID)
    {}

    virtual void outputXML(magnet::xml::XmlStream&) const = 0;
  
    shared_ptr<Property> _mass;
    shared_ptr<IDRange> range;
    std::string spName;
    unsigned int ID;
  };

  inline magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream& XML, const Species& g)
  {
    g.outputXML(XML);
    return XML;
  }
}
