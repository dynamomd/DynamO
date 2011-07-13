/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
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
#include "../../dynamics/ranges/1range.hpp"
#include "../../dynamics/globals/neighbourList.hpp"
#include <magnet/cloneptr.hpp>

namespace magnet { namespace xml { class Node; } }
namespace xml { class XmlStream; }

class Particle;

class CSCEntry: public dynamo::SimBase
{
public:
  CSCEntry(dynamo::SimData* const, const char *);
  
  virtual ~CSCEntry() {};

  virtual void initialise() {};
  
  friend magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream&, const CSCEntry&);

  static CSCEntry* getClass(const magnet::xml::Node&, dynamo::SimData* const);
 
  virtual void operator<<(const magnet::xml::Node&) = 0;

  bool isApplicable(const Particle& part) const;
  
  virtual void getParticleNeighbourhood(const Particle&, 
					const GNeighbourList::nbHoodFunc&) const {}

  virtual void getParticleLocalNeighbourhood(const Particle&, 
					     const GNeighbourList::nbHoodFunc&
					     ) const 
  {}

  virtual CSCEntry* Clone() const = 0;

protected:

  virtual void outputXML(magnet::xml::XmlStream&) const = 0;
  
  magnet::ClonePtr<CRange> range;
};
