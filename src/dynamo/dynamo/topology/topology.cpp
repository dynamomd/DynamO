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

#include <cstring>
#include <dynamo/particle.hpp>
#include <dynamo/ranges/IDRange.hpp>
#include <dynamo/ranges/IDRangeAll.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/topology/include.hpp>
#include <dynamo/topology/topology.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {

Topology::Topology(dynamo::Simulation *tmp, size_t nID)
    : SimBase_const(tmp, "Species"), ID(nID) {}

magnet::xml::XmlStream &operator<<(magnet::xml::XmlStream &XML,
                                   const Topology &g) {
  g.outputXML(XML);
  return XML;
}

void Topology::operator<<(const magnet::xml::Node &XML) {
  _name = XML.getAttribute("Name");

  if (!XML.hasNode("Molecule"))
    M_throw() << "Cannot load a Topology which has no molecules!";
}

shared_ptr<Topology> Topology::getClass(const magnet::xml::Node &XML,
                                        dynamo::Simulation *Sim, size_t ID) {
  if (!XML.getAttribute("Type").getValue().compare("Chain"))
    return shared_ptr<Topology>(new TChain(XML, Sim, ID));
  else if (!XML.getAttribute("Type").getValue().compare("PRIME"))
    return shared_ptr<Topology>(new TPRIME(XML, Sim, ID));
  else
    M_throw() << XML.getAttribute("Type").getValue()
              << ", Unknown type of Topology encountered";
}
} // namespace dynamo
