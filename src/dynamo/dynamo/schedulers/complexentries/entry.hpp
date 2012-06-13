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
#include <dynamo/ranges/1range.hpp>
#include <dynamo/globals/neighbourList.hpp>

namespace magnet { namespace xml { class Node; } }
namespace xml { class XmlStream; }

namespace dynamo {
  class Particle;

  class SCEntry: public dynamo::SimBase
  {
  public:
    SCEntry(dynamo::SimData* const, const char *);
  
    virtual ~SCEntry() {};

    virtual void initialise() {};
  
    friend magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream&, const SCEntry&);

    static SCEntry* getClass(const magnet::xml::Node&, dynamo::SimData* const);
 
    virtual void operator<<(const magnet::xml::Node&) = 0;

    bool isApplicable(const Particle& part) const;
  
    virtual void getParticleNeighbourhood(const Vector&, 
					  const GNeighbourList::nbHoodFunc2&) const = 0;

    virtual void getParticleNeighbourhood(const Particle&, 
					  const GNeighbourList::nbHoodFunc&) const = 0;

    virtual void getLocalNeighbourhood(const Particle&, 
				       const GNeighbourList::nbHoodFunc&
				       ) const = 0;
  protected:

    virtual void outputXML(magnet::xml::XmlStream&) const = 0;
  
    shared_ptr<Range> range;
  };
}
