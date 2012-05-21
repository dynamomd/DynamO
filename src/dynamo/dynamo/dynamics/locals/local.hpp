/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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

#include <dynamo/base.hpp>
#include <dynamo/dynamics/ranges/1range.hpp>
#include <magnet/math/vector.hpp>
#include <string>

namespace magnet { namespace xml { class Node; } }
namespace xml { class XmlStream; }
namespace dynamo {
  class IntEvent;
  class NEventData;
  class LocalEvent;

  /*! \brief Represents 1-particle event sources which are Local in space.
   *
   * The purpose of this specialized class is to allow 1-particle
   * events, which are localized in space, to be inserted into a
   * neighbor list for efficiency.
   *
   * To do this, the Local class provides the isInCell method, used by a
   * GNeighbourList to check if this Local is in a certain cell.
   */
  class Local: public dynamo::SimBase
  {
  public:
    Local(dynamo::SimData*, const char *);

    Local(Range*, dynamo::SimData*, const char *);
  
    virtual ~Local() {}

    bool isInteraction(const Particle&) const;

    virtual LocalEvent getEvent(const Particle&) const = 0;

    virtual void runEvent(Particle&, const LocalEvent&) const = 0;
  
    virtual bool isInCell(const Vector &, const Vector &) const = 0;

    virtual void initialise(size_t) = 0;

    friend magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream&, const Local&);

    static shared_ptr<Local> getClass(const magnet::xml::Node&, dynamo::SimData*);

    virtual void operator<<(const magnet::xml::Node&) = 0;

    void setName(const std::string& tmp) { localName = tmp; }

    const std::string& getName() const { return localName; }

    inline const size_t& getID() const { return ID; }

    virtual void checkOverlaps(const Particle&) const  {}

  protected:
    virtual void outputXML(magnet::xml::XmlStream&) const = 0;

    shared_ptr<Range> range;  
    std::string localName;
    size_t ID;
  };
}
