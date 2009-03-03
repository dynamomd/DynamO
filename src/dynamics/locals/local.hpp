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

#ifndef CLocal_HPP
#define CLocal_HPP

#include <string>
#include "../../datatypes/pluginpointer.hpp"
#include "../../base/is_base.hpp"
#include "../ranges/1range.hpp"
#include "../../datatypes/vector.hpp"

class XMLNode;
namespace xmlw
{
  class XmlStream;
}
class CIntEvent;
class CNParticleData;
class CLocalEvent;

class CLocal: public DYNAMO::SimBase
{
public:
  CLocal(DYNAMO::SimData*, const char *);

  CLocal(CRange*, DYNAMO::SimData*, const char *);
  
  virtual ~CLocal() {}

  bool isInteraction(const CParticle&) const;

  virtual CLocal* Clone() const = 0; //{ return new COPBlank(*this); };

  virtual CLocalEvent getEvent(const CParticle&) const = 0;

  virtual void runEvent(const CParticle&, const CLocalEvent&) const = 0;
  
  virtual bool isInCell(const CVector<>&, const CVector<>&) const = 0;

  virtual void initialise(size_t) = 0;

  friend xmlw::XmlStream& operator<<(xmlw::XmlStream&, const CLocal&);

  static CLocal* getClass(const XMLNode&, DYNAMO::SimData*);

  virtual void operator<<(const XMLNode&) = 0;

  void setName(const std::string& tmp) { localName = tmp; }

  const std::string& getName() const { return localName; }

  inline const size_t& getID() const { return ID; }
  
protected:
  virtual void outputXML(xmlw::XmlStream&) const = 0;

  smrtPlugPtr<CRange> range;  
  std::string localName;
  size_t ID;
};

#endif
