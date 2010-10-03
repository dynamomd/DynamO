/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include "topology.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../ranges/1range.hpp"
#include "../ranges/1RAll.hpp"
#include "../../simulation/particle.hpp"
#include "../../base/is_simdata.hpp"
#include "include.hpp"

Topology::Topology(DYNAMO::SimData* tmp, size_t nID):
  SimBase_const(tmp,"Species", IC_blue),
  ID(nID)
{ }

xml::XmlStream& operator<<(xml::XmlStream& XML, const Topology& g)
{
  g.outputXML(XML);
  return XML;
}

void 
Topology::operator<<(const XMLNode& XML)
{
    try {
      spName = XML.getAttribute("Name");
    } 
    catch (boost::bad_lexical_cast &)
      {
	M_throw() << "Failed a lexical cast in CTopology";
      }
    
    for (int i = 0; i < XML.nChildNode(); i++)
      ranges.push_back(magnet::ClonePtr<CRange>(CRange::loadClass(XML.getChildNode(i), Sim)));
}

void
Topology::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Name") << spName;
  
  BOOST_FOREACH(const magnet::ClonePtr<CRange>& plugPtr, ranges)
    XML << xml::tag("Molecule") << plugPtr
	<< xml::endtag("Molecule");
}


Topology* 
Topology::loadClass(const XMLNode& XML, DYNAMO::SimData* Sim, size_t ID)
{
  if (!strcmp(XML.getAttribute("Type"),"Chain"))
    return new CTChain(XML, Sim, ID);
  else 
    M_throw() << XML.getAttribute("Type")
	      << ", Unknown type of Topology encountered";
}
