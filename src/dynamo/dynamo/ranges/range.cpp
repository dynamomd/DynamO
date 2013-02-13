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
    if (!strcmp(XML.getAttribute("Type"),"All"))
      return new IDRangeAll(XML, Sim);
    else if (!strcmp(XML.getAttribute("Type"),"None"))
      return new IDRangeNone(XML);
    else if (!strcmp(XML.getAttribute("Type"),"Ranged"))
      return new IDRangeRange(XML);
    else if (!strcmp(XML.getAttribute("Type"),"List"))
      return new IDRangeList(XML);
    else if (!strcmp(XML.getAttribute("Type"),"Union"))
      return new IDRangeUnion(XML, Sim);
    else
      M_throw() << "Unknown type of IDRange encountered (" << XML.getAttribute("Type") << ")";
  }

  IDPairRange*
  IDPairRange::getClass(const magnet::xml::Node& XML, const dynamo::Simulation* Sim)
  {
    if (!strcmp(XML.getAttribute("Type"),"Pair"))
      return new IDPairRangePair(XML, Sim);
    else if (!strcmp(XML.getAttribute("Type"),"List"))
      return new IDPairRangeList(XML);
    else if (!strcmp(XML.getAttribute("Type"),"Single"))
      return new IDPairRangeSingle(XML,Sim);
    else if (!strcmp(XML.getAttribute("Type"),"Union"))
      return new IDPairRangeUnion(XML,Sim);
    else if (!strcmp(XML.getAttribute("Type"),"Chains"))
      return new IDPairRangeChains(XML,Sim);              
    else if (!strcmp(XML.getAttribute("Type"),"ChainGroups"))
      return new IDPairRangeChainGroups(XML,Sim);
    else if (!strcmp(XML.getAttribute("Type"),"ChainEnds"))
      return new IDPairRangeChainEnds(XML,Sim);
    else if (!strcmp(XML.getAttribute("Type"),"IntraChains"))
      return new IDPairRangeIntraChains(XML,Sim);              
    else if (!strcmp(XML.getAttribute("Type"),"Rings"))
      return new IDPairRangeRings(XML,Sim);
    else if (!strcmp(XML.getAttribute("Type"),"All"))
      return new IDPairRangeAll(XML,Sim);
    else if (!strcmp(XML.getAttribute("Type"),"None"))
      return new IDPairRangeNone(XML,Sim);
    else 
      M_throw() << "Unknown type of IDPairRange encountered (" << XML.getAttribute("Type") << ")";
  }
}
