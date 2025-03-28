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
#include <dynamo/interactions/include.hpp>
#include <dynamo/interactions/interaction.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/species/species.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
Interaction::Interaction(dynamo::Simulation *tmp, IDPairRange *nR)
    : SimBase(tmp, "Interaction"), range(nR) {}

void Interaction::operator<<(const magnet::xml::Node &XML) {
  range = shared_ptr<IDPairRange>(
      IDPairRange::getClass(XML.getNode("IDPairRange"), Sim));
  intName = XML.getAttribute("Name");
}

bool Interaction::isInteraction(const Event &coll) const {
  return isInteraction(Sim->particles[coll._particle1ID],
                       Sim->particles[coll._particle2ID]);
}

magnet::xml::XmlStream &operator<<(magnet::xml::XmlStream &XML,
                                   const Interaction &g) {
  g.outputXML(XML);
  return XML;
}

shared_ptr<IDPairRange> &Interaction::getRange() { return range; }

const shared_ptr<IDPairRange> &Interaction::getRange() const { return range; }

shared_ptr<Interaction> Interaction::getClass(const magnet::xml::Node &XML,
                                              dynamo::Simulation *Sim) {
  if (!XML.getAttribute("Type").getValue().compare("HardSphere"))
    return shared_ptr<Interaction>(new IHardSphere(XML, Sim));
  else if (!XML.getAttribute("Type").getValue().compare("SquareWell"))
    return shared_ptr<Interaction>(new ISquareWell(XML, Sim));
  else if (!XML.getAttribute("Type").getValue().compare("PRIME"))
    return shared_ptr<Interaction>(new IPRIME(XML, Sim));
  else if (!XML.getAttribute("Type").getValue().compare("ThinThread"))
    return shared_ptr<Interaction>(new IThinThread(XML, Sim));
  else if (!XML.getAttribute("Type").getValue().compare("SquareWellSeq"))
    return shared_ptr<Interaction>(new ISWSequence(XML, Sim));
  else if (!XML.getAttribute("Type").getValue().compare("SquareBond"))
    return shared_ptr<Interaction>(new ISquareBond(XML, Sim));
  else if (!XML.getAttribute("Type").getValue().compare("Null"))
    return shared_ptr<Interaction>(new INull(XML, Sim));
  else if (!XML.getAttribute("Type").getValue().compare("Lines"))
    return shared_ptr<Interaction>(new ILines(XML, Sim));
  else if (!XML.getAttribute("Type").getValue().compare("DSMC"))
    return shared_ptr<Interaction>(new IDSMC(XML, Sim));
  else if (!XML.getAttribute("Type").getValue().compare("Dumbbells"))
    return shared_ptr<Interaction>(new IDumbbells(XML, Sim));
  else if (!XML.getAttribute("Type").getValue().compare("ParallelCubes"))
    return shared_ptr<Interaction>(new IParallelCubes(XML, Sim));
  else if (!XML.getAttribute("Type").getValue().compare("Stepped"))
    return shared_ptr<Interaction>(new IStepped(XML, Sim));
  else
    M_throw() << XML.getAttribute("Type").getValue()
              << ", Unknown type of interaction encountered";
}
} // namespace dynamo
