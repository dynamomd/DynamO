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
  IDRange* 
  IDRange::getClass(const magnet::xml::Node& XML, const dynamo::Simulation * Sim)
  {
    if (!strcmp(XML.getAttribute("Range"),"All"))
      return new IDRangeAll(XML, Sim);
    else if (!strcmp(XML.getAttribute("Range"),"None"))
      return new IDRangeNone(XML);
    else if (!strcmp(XML.getAttribute("Range"),"Single"))
      return new IDRangeSingle(XML);
    else if (!strcmp(XML.getAttribute("Range"),"Ranged"))
      return new IDRangeRange(XML);
    else if (!strcmp(XML.getAttribute("Range"),"List"))
      return new IDRangeList(XML);
    else 
      M_throw() << XML.getAttribute("Range")
		<< ", Unknown type of IDRange encountered";
  }

  IDPairRange*
  IDPairRange::getClass(const magnet::xml::Node& XML, const dynamo::Simulation* Sim)
  {
    if (!strcmp(XML.getAttribute("Range"),"Pair"))
      return new IDPairRangePair(XML, Sim);
    else if (!strcmp(XML.getAttribute("Range"),"List"))
      return new IDPairRangeList(XML);
    else if (!strcmp(XML.getAttribute("Range"),"2Single"))
      return new IDPairRangeSingle(XML,Sim);
    else if (!strcmp(XML.getAttribute("Range"),"RangeList"))
      return new IDPairRangeRangeList(XML,Sim);
    else if (!strcmp(XML.getAttribute("Range"),"Chains"))
      return new IDPairRangeChains(XML,Sim);              
    else if (!strcmp(XML.getAttribute("Range"),"ChainGroups"))
      return new IDPairRangeChainGroups(XML,Sim);
    else if (!strcmp(XML.getAttribute("Range"),"ChainEnds"))
      return new IDPairRangeChainEnds(XML,Sim);
    else if (!strcmp(XML.getAttribute("Range"),"IntraChains"))
      return new IDPairRangeIntraChains(XML,Sim);              
    else if (!strcmp(XML.getAttribute("Range"),"Rings"))
      return new IDPairRangeRings(XML,Sim);
    else if (!strcmp(XML.getAttribute("Range"),"2All"))
      return new IDPairRangeAll(XML,Sim);
    else if (!strcmp(XML.getAttribute("Range"),"2None"))
      return new IDPairRangeNone(XML,Sim);
    else 
      M_throw() << XML.getAttribute("Range")
		<< ", Unknown type of IDPairRange encountered";
  }
}
