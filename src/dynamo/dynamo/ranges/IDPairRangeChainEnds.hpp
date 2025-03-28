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
#include <dynamo/particle.hpp>
#include <dynamo/ranges/IDPairRange.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
class IDPairRangeChainEnds : public IDPairRange {
public:
  IDPairRangeChainEnds(const magnet::xml::Node &XML, const dynamo::Simulation *)
      : rangeStart(0), rangeEnd(0), interval(0) {
    rangeStart = XML.getAttribute("Start").as<size_t>();
    rangeEnd = XML.getAttribute("End").as<size_t>();
    interval = XML.getAttribute("Interval").as<size_t>();

    // Guarrantee that they are ordered
    if (rangeStart > rangeEnd)
      std::swap(rangeStart, rangeEnd);

    if ((rangeEnd - rangeStart + 1) % interval)
      M_throw() << "Length of range does not split into an integer"
                << " number of intervals";
  }

  IDPairRangeChainEnds(size_t r1, size_t r2, size_t l)
      : rangeStart(r1), rangeEnd(r2), interval(l) {
    // Guarrantee that they are ordered
    if (rangeStart > rangeEnd)
      std::swap(rangeStart, rangeEnd);

    if ((rangeEnd - rangeStart + 1) % interval)
      M_throw() << "Length of range does not split into an integer"
                << " number of intervals";
  }

  virtual bool isInRange(const Particle &p1, const Particle &p2) const {
    if (p1.getID() > p2.getID())
      return ((p1.getID() <= rangeEnd) && (p2.getID() >= rangeStart) &&
              !((p2.getID() - rangeStart) % interval) &&
              (p1.getID() - p2.getID() == interval - 1));
    else
      return ((p2.getID() <= rangeEnd) && (p1.getID() >= rangeStart) &&
              !((p1.getID() - rangeStart) % interval) &&
              (p2.getID() - p1.getID() == interval - 1));
  }

  virtual bool isInRange(const Particle &p1) const {
    return (p1.getID() <= rangeEnd) &&
           (p1.getID() >= rangeStart) // Check if the particle is in the range
           && (!((p1.getID() - rangeStart) % interval) // is it the start?
               || !((p1.getID() - rangeStart + 1) % interval)); // Or the end?
  }

protected:
  virtual void outputXML(magnet::xml::XmlStream &XML) const {
    XML << magnet::xml::attr("Type") << "ChainEnds"
        << magnet::xml::attr("Start") << rangeStart << magnet::xml::attr("End")
        << rangeEnd << magnet::xml::attr("Interval") << interval;
  }

  size_t rangeStart;
  size_t rangeEnd;
  size_t interval;
};
} // namespace dynamo
