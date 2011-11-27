/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
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

#include <dynamo/dynamics/ranges/1RRange.hpp>
#include <dynamo/simulation/particle.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
  RRange::RRange(const magnet::xml::Node& XML) 
  { operator<<(XML); }

  void 
  RRange::operator<<(const magnet::xml::Node& XML)
  {
    if (strcmp(XML.getAttribute("Range"), "Ranged"))
      M_throw() << "Attempting to load RRange from non range";

    try {
      startID = XML.getAttribute("Start").as<unsigned long>();
      endID = XML.getAttribute("End").as<unsigned long>();
    }
    catch (boost::bad_lexical_cast &)
      { M_throw() << "Failed a lexical cast in RRange"; }
  }

  void 
  RRange::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Range") << "Ranged"
	<< magnet::xml::attr("Start") << startID
	<< magnet::xml::attr("End") << endID;
  }
}
