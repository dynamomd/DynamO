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

#include <dynamo/dynamics/ranges/2RSingle.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
  C2RSingle::C2RSingle(const magnet::xml::Node& XML, const dynamo::SimData* Sim)
  { 
    if (strcmp(XML.getAttribute("Range"),"2Single"))
      M_throw() << "Attempting to load a 2Single from a non pair";
  
    range = std::tr1::shared_ptr<CRange>(CRange::getClass(XML.getNode("SingleRange"), Sim));
  }

  bool 
  C2RSingle::isInRange(const Particle&p1, const Particle&p2) const
  {
    return (range->isInRange(p1) && range->isInRange(p2));
  }

  void 
  C2RSingle::operator<<(const magnet::xml::Node&)
  {
    M_throw() << "Due to problems with C2RSingle operator<< cannot work for this class";
  }

  void 
  C2RSingle::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Range") << "2Single" 
	<< magnet::xml::tag("SingleRange")
	<< range
	<< magnet::xml::endtag("SingleRange");
  }
}
