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
#include <dynamo/species/include.hpp>
#include <dynamo/units/units.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
void SpPoint::operator<<(const magnet::xml::Node &XML) {
  range = shared_ptr<IDRange>(IDRange::getClass(XML.getNode("IDRange"), Sim));
  _mass = Sim->_properties.getProperty(XML.getAttribute("Mass"),
                                       Property::Units::Mass());
  spName = XML.getAttribute("Name");
}

void SpPoint::outputXML(magnet::xml::XmlStream &XML) const {
  XML << magnet::xml::attr("Mass") << _mass->getName()
      << magnet::xml::attr("Name") << spName << magnet::xml::attr("Type")
      << "Point" << *range;
}

double SpPoint::getParticleKineticEnergy(size_t ID) const {
  const Particle &part = Sim->particles[ID];

#ifdef DYNAMO_DEBUG
  if (!isSpecies(part))
    M_throw() << "Getting the energy of a particle which does not belong to "
                 "this Species!";
#endif

  const double mass = getMass(ID);
  if (std::isinf(mass))
    return 0;

  if (std::dynamic_pointer_cast<BCLeesEdwards>(Sim->BCs))
    return 0.5 *
           static_cast<const BCLeesEdwards &>(*Sim->BCs)
               .getPeculiarVelocity(part)
               .nrm2() *
           mass;
  else
    return 0.5 * part.getVelocity().nrm2() * mass;
}
} // namespace dynamo
