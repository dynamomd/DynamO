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
#include <magnet/math/vector.hpp>

namespace magnet { namespace xml { class Node; } }
namespace xml { class XmlStream; }
class IntEvent;
class NEventData;
class LocalEvent;

class Local: public DYNAMO::SimBase
{
public:
  Local(DYNAMO::SimData*, const char *);

  Local(CRange*, DYNAMO::SimData*, const char *);
  
  virtual ~Local() {}

  bool isInteraction(const Particle&) const;

  virtual Local* Clone() const = 0; //{ return new OPBlank(*this); };

  virtual LocalEvent getEvent(const Particle&) const = 0;

  virtual void runEvent(const Particle&, const LocalEvent&) const = 0;
  
  virtual bool isInCell(const Vector &, const Vector &) const = 0;

  virtual void initialise(size_t) = 0;

  friend xml::XmlStream& operator<<(xml::XmlStream&, const Local&);

  static Local* getClass(const magnet::xml::Node&, DYNAMO::SimData*);

  virtual void operator<<(const magnet::xml::Node&) = 0;

  void setName(const std::string& tmp) { localName = tmp; }

  const std::string& getName() const { return localName; }

  inline const size_t& getID() const { return ID; }

  virtual void write_povray_info(std::ostream&) const {}

  virtual void checkOverlaps(const Particle&) const  {}

protected:
  virtual void outputXML(xml::XmlStream&) const = 0;

  magnet::ClonePtr<CRange> range;  
  std::string localName;
  size_t ID;
};
