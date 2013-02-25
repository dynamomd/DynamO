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

#include <dynamo/ranges/include.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
  magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream& XML,
				     const IDRange& range)
  { 
    XML << magnet::xml::tag("IDRange");
    range.outputXML(XML);
    XML << magnet::xml::endtag("IDRange");
    return XML;
  }

  magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream& XML,
				     const IDPairRange& range)
  { 
    XML << magnet::xml::tag("IDPairRange");
    range.outputXML(XML); 
    XML << magnet::xml::endtag("IDPairRange");
    return XML; 
  }

  IDRange* 
  IDRange::getClass(const magnet::xml::Node& XML, const dynamo::Simulation * Sim)
  {
    if (!XML.getAttribute("Type").getValue().compare("All"))
      return new IDRangeAll(XML, Sim);
    else if (!XML.getAttribute("Type").getValue().compare("None"))
      return new IDRangeNone(XML);
    else if (!XML.getAttribute("Type").getValue().compare("Ranged"))
      return new IDRangeRange(XML);
    else if (!XML.getAttribute("Type").getValue().compare("List"))
      return new IDRangeList(XML);
    else if (!XML.getAttribute("Type").getValue().compare("Union"))
      return new IDRangeUnion(XML, Sim);
    else
      M_throw() << "Unknown type of IDRange encountered (" << XML.getAttribute("Type").getValue() << ")";
  }

  IDPairRange*
  IDPairRange::getClass(const magnet::xml::Node& XML, const dynamo::Simulation* Sim)
  {
    if (!XML.getAttribute("Type").getValue().compare("Pair"))
      return new IDPairRangePair(XML, Sim);
    else if (!XML.getAttribute("Type").getValue().compare("List"))
      return new IDPairRangeList(XML);
    else if (!XML.getAttribute("Type").getValue().compare("Single"))
      return new IDPairRangeSingle(XML,Sim);
    else if (!XML.getAttribute("Type").getValue().compare("Union"))
      return new IDPairRangeUnion(XML,Sim);
    else if (!XML.getAttribute("Type").getValue().compare("Chains"))
      return new IDPairRangeChains(XML,Sim);              
    else if (!XML.getAttribute("Type").getValue().compare("ChainGroups"))
      return new IDPairRangeChainGroups(XML,Sim);
    else if (!XML.getAttribute("Type").getValue().compare("ChainEnds"))
      return new IDPairRangeChainEnds(XML,Sim);
    else if (!XML.getAttribute("Type").getValue().compare("IntraChains"))
      return new IDPairRangeIntraChains(XML,Sim);              
    else if (!XML.getAttribute("Type").getValue().compare("Rings"))
      return new IDPairRangeRings(XML,Sim);
    else if (!XML.getAttribute("Type").getValue().compare("All"))
      return new IDPairRangeAll(XML,Sim);
    else if (!XML.getAttribute("Type").getValue().compare("None"))
      return new IDPairRangeNone(XML,Sim);
    else 
      M_throw() << "Unknown type of IDPairRange encountered (" << XML.getAttribute("Type").getValue() << ")";
  }
}
