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
#include <dynamo/species/include.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
Species::~Species() {}

shared_ptr<Species> Species::getClass(const magnet::xml::Node &XML,
                                      dynamo::Simulation *tmp, size_t nID) {
  if (!XML.getAttribute("Type").getValue().compare("Point"))
    return shared_ptr<Species>(new SpPoint(XML, tmp, nID));
  else if (!XML.getAttribute("Type").getValue().compare("SphericalTop") ||
           !XML.getAttribute("Type").getValue().compare("Lines"))
    return shared_ptr<Species>(new SpSphericalTop(XML, tmp, nID));
  else if (!XML.getAttribute("Type").getValue().compare("FixedCollider"))
    return shared_ptr<Species>(new SpFixedCollider(XML, tmp, nID));
  else
    M_throw() << XML.getAttribute("Type").getValue()
              << ", Unknown type of species encountered";
}

double Species::getParticleKineticEnergy(const Particle &p) const {
  return getParticleKineticEnergy(p.getID());
}
} // namespace dynamo
