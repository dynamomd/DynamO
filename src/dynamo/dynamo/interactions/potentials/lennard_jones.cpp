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

#include <dynamo/interactions/potentials/lennard_jones.hpp>
#include <dynamo/simulation.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
  void 
  PotentialLennardJones::outputXML(magnet::xml::XmlStream& XML) const {
    XML << magnet::xml::attr("Type") << "LennardJones";
  }

  void 
  PotentialLennardJones::operator<<(const magnet::xml::Node& XML) {
  }

  void 
  PotentialLennardJones::calculateNextStep() const {
    M_throw() << "Not implemented";
  }
}
