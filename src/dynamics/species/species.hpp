/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef CSPECIES_H
#define CSPECIES_H

#include "../../datatypes/pluginpointer.hpp"
#include "../../base/is_base.hpp"
#include "../ranges/1range.hpp"
#include "../units/units.hpp"
#include <string>
#include <list>

class XMLNode;
namespace xmlw
{
  class XmlStream;
}
class CParticle;
class CInteraction;

class CSpecHasInertia {};

class CSpecies:public DYNAMO::SimBase_const
{
public:  
  CSpecies(DYNAMO::SimData*, CRange*, Iflt nMass, std::string nName, 
	   unsigned int ID, std::string nIName="Bulk");
  
  CSpecies(const XMLNode&, DYNAMO::SimData*, unsigned int ID);

  virtual ~CSpecies() {}

  bool isSpecies(const CParticle &) const;
  
  const Iflt& getMass() const { return mass; }
  
  unsigned long getCount() const;
  
  unsigned int getID() const { return ID; }
  
  virtual void operator<<(const XMLNode&);

  void initialise();

  friend xmlw::XmlStream& operator<<(xmlw::XmlStream&, const CSpecies&);
  
  const std::string& getName() const { return spName; }

  const std::string& getIntName() const { return intName; }

  const CInteraction* getIntPtr() const;

  void setIntPtr(CInteraction*);

  const smrtPlugPtr<CRange>& getRange() const { return range; }

  virtual CSpecies* Clone() const { return new CSpecies(*this); }

  virtual Iflt getScalarMomentOfInertia() const { return 0; }

  static CSpecies* getClass(const XMLNode&, DYNAMO::SimData*, unsigned int);

protected:
  CSpecies(DYNAMO::SimData*, const char* name, const char* color, 
	   CRange*, Iflt nMass, std::string nName, 
	   unsigned int ID, std::string nIName="Bulk");
  

  virtual void outputXML(xmlw::XmlStream&) const;
  
  Iflt mass;

  smrtPlugPtr<CRange> range;

  std::string spName;
  std::string intName;

  CInteraction* IntPtr;

  unsigned int ID;
};



#endif
