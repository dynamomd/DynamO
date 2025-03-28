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
#include <dynamo/BC/include.hpp>
#include <dynamo/simulation.hpp>
#include <magnet/exception.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
magnet::xml::XmlStream &operator<<(magnet::xml::XmlStream &XML,
                                   const BoundaryCondition &g) {
  g.outputXML(XML);
  return XML;
}

namespace {
typedef shared_ptr<BoundaryCondition> retptr;
}

double BoundaryCondition::getDistance(const Particle &p1, const Particle &p2) {
  Vector rij = p1.getPosition() - p2.getPosition();
  applyBC(rij);
  return rij.nrm();
}

retptr BoundaryCondition::getClass(const magnet::xml::Node &XML,
                                   dynamo::Simulation *tmp) {
  if (!XML.getAttribute("Type").getValue().compare("None") ||
      !XML.getAttribute("Type").getValue().compare("Null"))
    return retptr(new BCNone(tmp));
  else if (!XML.getAttribute("Type").getValue().compare("PBC"))
    return retptr(new BCPeriodic(tmp));
  else if (!XML.getAttribute("Type").getValue().compare("NoXPBC"))
    return retptr(new BCPeriodicExceptX(tmp));
  else if (!XML.getAttribute("Type").getValue().compare("OnlyXPBC"))
    return retptr(new BCPeriodicXOnly(tmp));
  else if (!XML.getAttribute("Type").getValue().compare("LE"))
    return retptr(new BCLeesEdwards(XML, tmp));
  else
    M_throw() << XML.getAttribute("Type").getValue()
              << ", Unknown type of boundary encountered";
}
} // namespace dynamo
