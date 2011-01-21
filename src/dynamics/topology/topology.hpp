/*  DYNAMO:- Event driven molecular dynamics simulator 
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

#include <magnet/cloneptr.hpp>
#include "../../base/is_base.hpp"
#include "../ranges/1range.hpp"
#include <string>
#include <list>

class XMLNode;
namespace xml
{
  class XmlStream;
}
class Particle;
class Interaction;

class Topology:public DYNAMO::SimBase_const
{
public:  
  virtual ~Topology() {}

  bool isInStructure(const Particle &) const;
  
  const size_t& getID() const { return ID; }
  
  void operator<<(const XMLNode&);

  virtual void initialise() {}

  friend xml::XmlStream& operator<<(xml::XmlStream&, const Topology&);
  
  const std::string& getName() const
  { return spName; }
  
  static Topology* loadClass(const XMLNode& ,DYNAMO::SimData*, size_t);

  virtual Topology* Clone() const = 0; //{ return new CTopology(this); }

  inline void addMolecule(CRange* ptr)
  { ranges.push_back(magnet::ClonePtr<CRange>(ptr)); }

  inline const std::list<magnet::ClonePtr<CRange> >& getMolecules() const
  { return ranges; }

  inline size_t getMoleculeCount() const { return ranges.size(); }

protected:
  Topology(DYNAMO::SimData*, size_t ID);

  virtual void outputXML(xml::XmlStream&) const;
  
  std::list<magnet::ClonePtr<CRange> > ranges;
  
  std::string spName;
  
  size_t ID;
};
