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

#ifndef CSystem_HPP
#define CSystem_HPP

#include <string>
#include "../../base/is_base.hpp"
#include "../eventtypes.hpp"

class XMLNode;
namespace xmlw
{
  class XmlStream;
}
class IntEvent;
class GlobalEvent;
class NEventData;

class CSystem: public DYNAMO::SimBase
{
public:
  CSystem(DYNAMO::SimData*);
  
  virtual ~CSystem() {}
  
  virtual CSystem* Clone() const = 0; //{ return new OPBlank(*this); };

  inline void stream(const Iflt& ndt) { dt -= ndt; }

  virtual void runEvent() const = 0;

  virtual void initialise(size_t) = 0;

  virtual void operator<<(const XMLNode&) = 0;

  bool operator<(const IntEvent&) const;

  bool operator<(const GlobalEvent&) const;

  bool operator<(const CSystem&) const;
  
  Iflt getdt() const { return dt; }
  
  friend xmlw::XmlStream& operator<<(xmlw::XmlStream&, const CSystem&);
  
  static CSystem* getClass(const XMLNode&, DYNAMO::SimData*);
  
  void setName(const std::string& tmp) { sysName = tmp; }

  const std::string& getName() const { return sysName; }

  EEventType getType() const { return type; }

  virtual void changeSystem(DYNAMO::SimData* ptr) { Sim = ptr; }

  inline const size_t& getID() const { return ID; }

protected:
  virtual void outputXML(xmlw::XmlStream&) const = 0;

  std::string sysName;  
  mutable Iflt dt;
  EEventType type;
  size_t ID;
};

#endif
