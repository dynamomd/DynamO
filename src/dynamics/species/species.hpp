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
#include "../coilRenderObj.hpp"

namespace magnet { namespace xml { class Node; } }
namespace xml { class XmlStream; }
class Particle;
class Interaction;
class RenderObj;
namespace magnet { namespace GL { class CLGLState; } }

class Species:public DYNAMO::SimBase, public CoilRenderObj
{
public:
  virtual ~Species() {}

  inline bool isSpecies(const Particle& p1) const { return range->isInRange(p1); }  

  inline const double& getMass() const { return mass; }
  inline unsigned long getCount() const { return range->size(); }
  inline unsigned int getID() const { return ID; }
  inline const std::string& getName() const { return spName; }
  inline const std::string& getIntName() const { return intName; }
  inline const magnet::ClonePtr<CRange>& getRange() const { return range; }
  inline const Interaction* getIntPtr() const { return IntPtr; }
  inline void setIntPtr(Interaction* nPtr) { IntPtr = nPtr; }
  virtual double getScalarMomentOfInertia() const = 0;

  virtual void operator<<(const magnet::xml::Node&) = 0;

  virtual void initialise() = 0;

  friend xml::XmlStream& operator<<(xml::XmlStream&, const Species&);
  
  virtual Species* Clone() const = 0;

  static Species* getClass(const magnet::xml::Node&, DYNAMO::SimData*, size_t);

protected:
  Species(DYNAMO::SimData* tmp, std::string name, 
	  CRange* nr, double nMass, std::string nName, 
	  unsigned int nID, std::string nIName):
    SimBase(tmp,name, IC_blue),
    mass(nMass),range(nr),spName(nName),intName(nIName),IntPtr(NULL),
    ID(nID)
  {}

  virtual void outputXML(xml::XmlStream&) const = 0;
  
  double mass;
  magnet::ClonePtr<CRange> range;
  std::string spName;
  std::string intName;
  Interaction* IntPtr;
  unsigned int ID;
};

inline xml::XmlStream& operator<<(xml::XmlStream& XML, const Species& g)
{
  g.outputXML(XML);
  return XML;
}
