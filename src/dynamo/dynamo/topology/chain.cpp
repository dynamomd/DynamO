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

#include <dynamo/topology/chain.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
  TChain::TChain(const magnet::xml::Node& XML, dynamo::Simulation* Sim, unsigned int ID):
    Topology(Sim, ID) 
  {
    TChain::operator<<(XML);
  }

  void 
  TChain::operator<<(const magnet::xml::Node& XML) 
  {
    Topology::operator<<(XML);

    for (magnet::xml::Node node = XML.findNode("Molecule"); node.valid(); ++node)
      ranges.push_back(shared_ptr<IDRange>(IDRange::getClass(node.getNode("IDRange"), Sim)));
  
    size_t Clength = (*ranges.begin())->size();
    for (const shared_ptr<IDRange>& nRange : ranges)
      if (nRange->size() != Clength)
	M_throw() << "Size mismatch in loading one of the ranges in Chain topology \"" 
		  << _name << "\"";
  }

  TChain::TChain(dynamo::Simulation* Sim, unsigned int ID, std::string nName):
    Topology(Sim,ID)
  {
    _name = nName;
  }

  void 
  TChain::outputXML(magnet::xml::XmlStream& XML) const 
  {
    XML << magnet::xml::attr("Name") << _name;
    XML << magnet::xml::attr("Type") << "Chain";
  
    for (const shared_ptr<IDRange>& plugPtr : ranges)
      XML << magnet::xml::tag("Molecule") << plugPtr
	  << magnet::xml::endtag("Molecule");
  }
}
