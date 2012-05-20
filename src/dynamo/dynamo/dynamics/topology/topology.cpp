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

#include <dynamo/dynamics/topology/topology.hpp>
#include <dynamo/dynamics/ranges/1range.hpp>
#include <dynamo/dynamics/ranges/1RAll.hpp>
#include <dynamo/simulation/particle.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <dynamo/dynamics/topology/include.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>

namespace dynamo {

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
    
    if (!XML.hasNode("Molecule"))
      M_throw() << "Cannot load a Topology which has no molecules!";

    for (magnet::xml::Node node = XML.fastGetNode("Molecule"); node.valid(); ++node)
      ranges.push_back(shared_ptr<Range>(Range::getClass(node, Sim)));
  }

  void
  Topology::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Name") << spName;
  
    BOOST_FOREACH(const shared_ptr<Range>& plugPtr, ranges)
      XML << magnet::xml::tag("Molecule") << plugPtr
	  << magnet::xml::endtag("Molecule");
  }


  shared_ptr<Topology>
  Topology::getClass(const magnet::xml::Node& XML, dynamo::SimData* Sim, size_t ID)
  {
    if (!strcmp(XML.getAttribute("Type"),"Chain"))
      return shared_ptr<Topology>(new TChain(XML, Sim, ID));
    else 
      M_throw() << XML.getAttribute("Type")
		<< ", Unknown type of Topology encountered";
  }
}
