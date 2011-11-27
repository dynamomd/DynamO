/*  dynamo:- Event driven molecular dynamics simulator 
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

#include <dynamo/dynamics/ranges/include.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
  Range* 
  Range::getClass(const magnet::xml::Node& XML, const dynamo::SimData * Sim)
  {
    if (!strcmp(XML.getAttribute("Range"),"All"))
      return new RAll(XML, Sim);
    else if (!strcmp(XML.getAttribute("Range"),"None"))
      return new RNone(XML);
    else if (!strcmp(XML.getAttribute("Range"),"Single"))
      return new RSingle(XML);
    else if (!strcmp(XML.getAttribute("Range"),"Ranged"))
      return new RRange(XML);
    else if (!strcmp(XML.getAttribute("Range"),"List"))
      return new RList(XML);
    else 
      M_throw() << XML.getAttribute("Range")
		<< ", Unknown type of Range encountered";
  }

  magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream& XML, const Range&g)
  {
    g.outputXML(XML);
    return XML;
  }

  magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream& XML, 
				     const C2Range& g)
  {
    g.outputXML(XML);
    return XML;
  }

  C2Range*
  C2Range::getClass(const magnet::xml::Node& XML, const dynamo::SimData* Sim)
  {
    if (!strcmp(XML.getAttribute("Range"),"Pair"))
      return new C2RPair(XML, Sim);
    else if (!strcmp(XML.getAttribute("Range"),"List"))
      return new C2RList(XML);
    else if (!strcmp(XML.getAttribute("Range"),"2Single"))
      return new C2RSingle(XML,Sim);
    else if (!strcmp(XML.getAttribute("Range"),"RangeList"))
      return new C2RRangeList(XML,Sim);
    else if (!strcmp(XML.getAttribute("Range"),"Chain"))
      return new C2RChain(XML,Sim);
    else if (!strcmp(XML.getAttribute("Range"),"Chains"))
      return new C2RChains(XML,Sim);              
    else if (!strcmp(XML.getAttribute("Range"),"ChainGroups"))
      return new C2RChainGroups(XML,Sim);
    else if (!strcmp(XML.getAttribute("Range"),"ChainEnds"))
      return new C2RChainEnds(XML,Sim);
    else if (!strcmp(XML.getAttribute("Range"),"IntraChains"))
      return new C2RIntraChains(XML,Sim);              
    else if (!strcmp(XML.getAttribute("Range"),"Ring"))
      return new C2RRing(XML,Sim);
    else if (!strcmp(XML.getAttribute("Range"),"Rings"))
      return new C2RRings(XML,Sim);
    else if (!strcmp(XML.getAttribute("Range"),"2All"))
      return new C2RAll(XML,Sim);
    else if (!strcmp(XML.getAttribute("Range"),"2None"))
      return new C2RNone(XML,Sim);
    else 
      M_throw() << XML.getAttribute("Range")
		<< ", Unknown type of C2Range encountered";
  }
}
