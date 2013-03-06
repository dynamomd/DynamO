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

#include <dynamo/topology/topology.hpp>
#include <dynamo/ranges/IDRange.hpp>
#include <dynamo/ranges/IDRangeAll.hpp>
#include <dynamo/particle.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/topology/include.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>

namespace dynamo {

  Topology::Topology(dynamo::Simulation* tmp, size_t nID):
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
    spName = XML.getAttribute("Name");
    
    if (!XML.hasNode("Molecule"))
      M_throw() << "Cannot load a Topology which has no molecules!";

    for (magnet::xml::Node node = XML.fastGetNode("Molecule"); node.valid(); ++node)
      ranges.push_back(shared_ptr<IDRange>(IDRange::getClass(node.getNode("IDRange"), Sim)));
  }

  void
  Topology::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Name") << spName;
  
    BOOST_FOREACH(const shared_ptr<IDRange>& plugPtr, ranges)
      XML << magnet::xml::tag("Molecule") << plugPtr
	  << magnet::xml::endtag("Molecule");
  }


  shared_ptr<Topology>
  Topology::getClass(const magnet::xml::Node& XML, dynamo::Simulation* Sim, size_t ID)
  {
    if (!XML.getAttribute("Type").getValue().compare("Chain"))
      return shared_ptr<Topology>(new TChain(XML, Sim, ID));
    else 
      M_throw() << XML.getAttribute("Type").getValue()
		<< ", Unknown type of Topology encountered";
  }
}
