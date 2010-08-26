/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "../../datatypes/pluginpointer.hpp"
#include "../../base/is_base.hpp"
#include "../ranges/1range.hpp"
#include "../units/units.hpp"
#include <string>
#include <list>

class XMLNode;
namespace xml
{
  class XmlStream;
}
class Particle;
class Interaction;


class Species:public DYNAMO::SimBase_const
{
public:  
  Species(DYNAMO::SimData*, CRange*, Iflt nMass, std::string nName, 
	   unsigned int ID, std::string nIName="Bulk");
  
  Species(const XMLNode&, DYNAMO::SimData*, unsigned int ID);

  virtual ~Species() {}

  bool isSpecies(const Particle &) const;
  
  const Iflt& getMass() const { return mass; }
  
  unsigned long getCount() const;
  
  unsigned int getID() const { return ID; }
  
  virtual void operator<<(const XMLNode&);

  void initialise();

  friend xml::XmlStream& operator<<(xml::XmlStream&, const Species&);
  
  const std::string& getName() const { return spName; }

  const std::string& getIntName() const { return intName; }

  const Interaction* getIntPtr() const;

  void setIntPtr(Interaction*);

  const ClonePtr<CRange>& getRange() const { return range; }

  virtual Species* Clone() const { return new Species(*this); }

  virtual Iflt getScalarMomentOfInertia() const { return 0; }

  static Species* getClass(const XMLNode&, DYNAMO::SimData*, unsigned int);

protected:
  Species(DYNAMO::SimData*, const char* name, const char* color, 
	   CRange*, Iflt nMass, std::string nName, 
	   unsigned int ID, std::string nIName="Bulk");
  

  virtual void outputXML(xml::XmlStream&) const;
  
  Iflt mass;

  ClonePtr<CRange> range;

  std::string spName;
  std::string intName;

  Interaction* IntPtr;

  unsigned int ID;
};

class SpInertia: public Species
{
public:
  SpInertia(DYNAMO::SimData* sim, const char* name, const char* color, 
	       CRange* r, Iflt nMass, std::string nName, 
	       unsigned int ID, std::string nIName="Bulk"):
  Species(sim, name, color, r, nMass, nName, ID, nIName)
  {}

  SpInertia(const XMLNode& XML, DYNAMO::SimData* Sim, unsigned int ID):
    Species(XML,Sim,ID)
  {}

};
