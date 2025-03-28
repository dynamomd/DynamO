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
#include <boost/functional/hash.hpp>
#include <dynamo/particle.hpp>
#include <dynamo/ranges/IDPairRange.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>
#include <unordered_set>

namespace dynamo {
class IDPairRangeList : public IDPairRange {
  typedef std::pair<unsigned long, unsigned long> Key;
  typedef std::unordered_set<Key, boost::hash<Key>> Container;

public:
  IDPairRangeList(const magnet::xml::Node &XML) { operator<<(XML); }

  IDPairRangeList() {}

  virtual bool isInRange(const Particle &p1, const Particle &p2) const {
    return pairmap.find(Key(std::min(p1.getID(), p2.getID()),
                            std::max(p1.getID(), p2.getID()))) != pairmap.end();
  }

  virtual bool isInRange(const Particle &p1) const {
    // Look for the key by hand!
    for (const auto &element : pairmap)
      if ((element.first == p1.getID()) || (element.second == p1.getID()))
        return true;

    return false;
  }

  void addPair(unsigned long a, unsigned long b) {
    pairmap.insert(Key(std::min(a, b), std::max(a, b)));
  }

  const Container &getPairMap() const { return pairmap; }

  virtual void operator<<(const magnet::xml::Node &XML) {
    for (magnet::xml::Node node = XML.findNode("IDPair"); node.valid(); ++node)
      addPair(node.getAttribute("ID1").as<unsigned long>(),
              node.getAttribute("ID2").as<unsigned long>());
  }

protected:
  virtual void outputXML(magnet::xml::XmlStream &XML) const {
    XML << magnet::xml::attr("Type") << "List";
    for (const Key &key : pairmap)
      XML << magnet::xml::tag("IDPair") << magnet::xml::attr("ID1") << key.first
          << magnet::xml::attr("ID2") << key.second
          << magnet::xml::endtag("IDPair");
  }

  Container pairmap;
};
} // namespace dynamo
