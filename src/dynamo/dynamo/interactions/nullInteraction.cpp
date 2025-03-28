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
#include <dynamo/2particleEventData.hpp>
#include <dynamo/interactions/nullInteraction.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
INull::INull(dynamo::Simulation *tmp, IDPairRange *nR, std::string name)
    : Interaction(tmp, nR) {
  intName = name;
}

INull::INull(const magnet::xml::Node &XML, dynamo::Simulation *tmp)
    : Interaction(tmp, NULL) {
  operator<<(XML);
}

void INull::initialise(size_t nID) { ID = nID; }

void INull::operator<<(const magnet::xml::Node &XML) {
  Interaction::operator<<(XML);
  intName = XML.getAttribute("Name");
}

Event INull::getEvent(const Particle &p1, const Particle &p2) const {
  return Event(p1, std::numeric_limits<float>::infinity(), INTERACTION, NONE,
               ID, p2);
}

PairEventData INull::runEvent(Particle &, Particle &, Event) {
  M_throw() << "Null event trying to run a collision!";
}

void INull::outputXML(magnet::xml::XmlStream &XML) const {
  XML << magnet::xml::attr("Type") << "Null" << magnet::xml::attr("Name")
      << intName << range;
}
} // namespace dynamo
