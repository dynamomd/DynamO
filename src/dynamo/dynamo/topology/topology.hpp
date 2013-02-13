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
#include <string>
#include <list>

namespace magnet { namespace xml { class Node; } }
namespace xml { class XmlStream; }
namespace dynamo {
  class Particle;
  class Interaction;

  class Topology: public dynamo::SimBase_const
  {
  public:  
    virtual ~Topology() {}

    bool isInStructure(const Particle &) const;
  
    const size_t& getID() const { return ID; }
  
    virtual void operator<<(const magnet::xml::Node&);

    virtual void initialise() {}

    friend magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream&, const Topology&);
  
    const std::string& getName() const
    { return spName; }
  
    static shared_ptr<Topology> getClass(const magnet::xml::Node& ,dynamo::Simulation*, size_t);

    inline void addMolecule(IDRange* ptr)
    { ranges.push_back(shared_ptr<IDRange>(ptr)); }

    inline const std::list<shared_ptr<IDRange> >& getMolecules() const
    { return ranges; }

    inline size_t getMoleculeCount() const { return ranges.size(); }

  protected:
    Topology(dynamo::Simulation*, size_t ID);

    virtual void outputXML(magnet::xml::XmlStream&) const;
  
    std::list<shared_ptr<IDRange> > ranges;
  
    std::string spName;
  
    size_t ID;
  };
}
