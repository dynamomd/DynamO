/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef CTopology_H
#define CTopology_H

#include "../../datatypes/pluginpointer.hpp"
#include "../../base/is_base.hpp"
#include "../ranges/1range.hpp"
#include <string>
#include <list>

class XMLNode;
namespace xmlw
{
  class XmlStream;
}
class CParticle;
class CInteraction;

class CTopology:public DYNAMO::SimBase_const
{
public:  
  virtual ~CTopology() {}

  bool isInStructure(const CParticle &) const;
  
  const size_t& getID() const { return ID; }
  
  void operator<<(const XMLNode&);

  virtual void initialise() {}

  friend xmlw::XmlStream& operator<<(xmlw::XmlStream&, const CTopology&);
  
  const std::string& getName() const
  { return spName; }
  
  static CTopology* loadClass(const XMLNode& ,DYNAMO::SimData*, size_t);

  virtual CTopology* Clone() const = 0; //{ return new CTopology(this); }

  inline void addMolecule(CRange* ptr)
  { ranges.push_back(smrtPlugPtr<CRange>(ptr)); }

  inline const std::list<smrtPlugPtr<CRange> >& getMolecules() const
  { return ranges; }

  inline size_t getMoleculeCount() const { return ranges.size(); }

protected:
  CTopology(DYNAMO::SimData*, size_t ID);

  virtual void outputXML(xmlw::XmlStream&) const;
  
  std::list<smrtPlugPtr<CRange> > ranges;
  
  std::string spName;
  
  size_t ID;
};

#endif
