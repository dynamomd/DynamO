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

#include <dynamo/species/include.hpp>
#include <dynamo/ranges/IDRange.hpp>
#include <dynamo/ranges/IDRangeAll.hpp>
#include <dynamo/particle.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/units/units.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <cstring>

namespace dynamo {
  void 
  SpPoint::operator<<(const magnet::xml::Node& XML)
  {
    range = shared_ptr<IDRange>(IDRange::getClass(XML.getNode("IDRange"), Sim));
    _mass = Sim->_properties.getProperty(XML.getAttribute("Mass"), Property::Units::Mass());
    spName = XML.getAttribute("Name");
  }

  void 
  SpPoint::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Mass") 
	<< _mass->getName()
	<< magnet::xml::attr("Name") << spName
	<< magnet::xml::attr("Type") << "Point"
	<< *range;
  }
}
