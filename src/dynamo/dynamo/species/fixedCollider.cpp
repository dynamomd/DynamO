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

#include <dynamo/species/fixedCollider.hpp>
#include <boost/foreach.hpp>
#include <dynamo/simulation.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
  void 
  SpFixedCollider::initialise()
  {
    SpPoint::initialise();

    BOOST_FOREACH(size_t ID, *range)
      Sim->particles[ID].clearState(Particle::DYNAMIC);
  }

  void 
  SpFixedCollider::operator<<(const magnet::xml::Node& XML)
  {
    range = shared_ptr<IDRange>(IDRange::getClass(XML.getNode("IDRange"), Sim));
    spName = XML.getAttribute("Name");
    intName = XML.getAttribute("IntName");
  }

  void 
  SpFixedCollider::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Name") << spName
	<< magnet::xml::attr("IntName") << intName
	<< magnet::xml::attr("Type") << "FixedCollider";
  
    XML << *range;
  }
}
