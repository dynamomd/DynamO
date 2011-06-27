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

#include "topology.hpp"
#include "../ranges/1range.hpp"
#include "../ranges/1RAll.hpp"
#include "../../simulation/particle.hpp"
#include "../../base/is_simdata.hpp"
#include "include.hpp"
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>

Topology::Topology(dynamo::SimData* tmp, size_t nID):
  SimBase_const(tmp, "Species"),
  ID(nID)
{ }

magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream& XML, const Topology& g)
{
  g.outputXML(XML);
  return XML;
}

void 
Topology::operator<<(const magnet::xml::Node& XML)
{
    try { spName = XML.getAttribute("Name"); } 
    catch (boost::bad_lexical_cast &)
      {
	M_throw() << "Failed a lexical cast in CTopology";
      }
    
    for (magnet::xml::Node node = XML.getNode("Molecule"); node.valid(); ++node)
      ranges.push_back(magnet::ClonePtr<CRange>(CRange::getClass(node, Sim)));
}

void
Topology::outputXML(magnet::xml::XmlStream& XML) const
{
  XML << magnet::xml::attr("Name") << spName;
  
  BOOST_FOREACH(const magnet::ClonePtr<CRange>& plugPtr, ranges)
    XML << magnet::xml::tag("Molecule") << plugPtr
	<< magnet::xml::endtag("Molecule");
}


Topology* 
Topology::getClass(const magnet::xml::Node& XML, dynamo::SimData* Sim, size_t ID)
{
  if (!strcmp(XML.getAttribute("Type"),"Chain"))
    return new CTChain(XML, Sim, ID);
  else 
    M_throw() << XML.getAttribute("Type")
	      << ", Unknown type of Topology encountered";
}
