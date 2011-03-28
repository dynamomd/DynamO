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
#include "../../base/is_base.hpp"
#include "../eventtypes.hpp"
#include <string>

struct XMLNode;
namespace xml
{
  class XmlStream;
}
class IntEvent;
class GlobalEvent;
class NEventData;

class System: public DYNAMO::SimBase
{
public:
  System(DYNAMO::SimData*);
  
  virtual ~System() {}
  
  virtual System* Clone() const = 0; //{ return new OPBlank(*this); };

  inline void stream(const double& ndt) { dt -= ndt; }

  virtual void runEvent() const = 0;

  virtual void initialise(size_t) = 0;

  virtual void operator<<(const XMLNode&) = 0;

  bool operator<(const IntEvent&) const;

  bool operator<(const GlobalEvent&) const;

  bool operator<(const System&) const;
  
  double getdt() const { return dt; }
  
  friend xml::XmlStream& operator<<(xml::XmlStream&, const System&);
  
  static System* getClass(const XMLNode&, DYNAMO::SimData*);
  
  void setName(const std::string& tmp) { sysName = tmp; }

  const std::string& getName() const { return sysName; }

  EEventType getType() const { return type; }

  virtual void changeSystem(DYNAMO::SimData* ptr) { Sim = ptr; }

  inline const size_t& getID() const { return ID; }

protected:
  virtual void outputXML(xml::XmlStream&) const = 0;

  std::string sysName;  
  mutable double dt;
  EEventType type;
  size_t ID;
};
