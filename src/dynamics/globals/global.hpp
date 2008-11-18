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

#ifndef CGlobal_HPP
#define CGlobal_HPP

#include <string>
#include "../../datatypes/pluginpointer.hpp"
#include "../../base/is_base.hpp"
#include "../ranges/1range.hpp"

class XMLNode;
namespace xmlw
{
  class XmlStream;
}
class CIntEvent;
class CNParticleData;
class CGlobEvent;

class CGlobal: public DYNAMO::SimBase
{
public:
  CGlobal(DYNAMO::SimData*, const char *);

  CGlobal(CRange*, DYNAMO::SimData*, const char *);
  
  virtual ~CGlobal() {}

  bool isInteraction(const CParticle &) const;

  virtual CGlobal* Clone() const = 0; //{ return new COPBlank(*this); };

  virtual CGlobEvent getEvent(const CParticle &) const = 0;

  virtual void runEvent(const CParticle&) const = 0;

  virtual void initialise(size_t) = 0;

  friend xmlw::XmlStream& operator<<(xmlw::XmlStream&, const CGlobal&);

  static CGlobal* getClass(const XMLNode&, DYNAMO::SimData*);

  virtual void operator<<(const XMLNode&) = 0;

  void setName(const std::string& tmp) { globName = tmp; }

  const std::string& getName() const { return globName; }

  inline const size_t& getID() const { return ID; }
  
protected:
  virtual void outputXML(xmlw::XmlStream&) const = 0;

  smrtPlugPtr<CRange> range;  
  std::string globName;
  size_t ID;
};

#endif
