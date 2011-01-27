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

#include <string>
#include <magnet/cloneptr.hpp>
#include "../../base/is_base.hpp"
#include "../ranges/1range.hpp"

class XMLNode;
namespace xml
{
  class XmlStream;
}
class IntEvent;
class NEventData;
class GlobalEvent;

class Global: public DYNAMO::SimBase
{
public:
  Global(DYNAMO::SimData*, const char *);

  Global(CRange*, DYNAMO::SimData*, const char *);
  
  virtual ~Global() {}

  bool isInteraction(const Particle &) const;

  virtual Global* Clone() const = 0; //{ return new OPBlank(*this); };

  virtual GlobalEvent getEvent(const Particle &) const = 0;

  virtual void runEvent(const Particle&, const double) const = 0;

  virtual void initialise(size_t) = 0;

  friend xml::XmlStream& operator<<(xml::XmlStream&, const Global&);

  static Global* getClass(const XMLNode&, DYNAMO::SimData*);

  virtual void operator<<(const XMLNode&) = 0;

  void setName(const std::string& tmp) { globName = tmp; }

  const std::string& getName() const { return globName; }

  inline const size_t& getID() const { return ID; }
  
protected:
  virtual void outputXML(xml::XmlStream&) const = 0;

  magnet::ClonePtr<CRange> range;  
  std::string globName;
  size_t ID;
};
