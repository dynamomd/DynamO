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

#include <dynamo/dynamics/ranges/2RPair.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
  C2RPair::C2RPair(const magnet::xml::Node& XML, const dynamo::SimData* Sim)
  { 
    if (strcmp(XML.getAttribute("Range"), "Pair"))
      M_throw() << "Attempting to load a pair from a non pair";
  
    range1 = shared_ptr<Range>(Range::getClass(XML.getNode("Range1"), Sim));
    range2 = shared_ptr<Range>(Range::getClass(XML.getNode("Range2"), Sim));
  }

  bool 
  C2RPair::isInRange(const Particle&p1, const Particle&p2) const
  {
    if ((range1->isInRange(p1) && range2->isInRange(p2))
	|| (range1->isInRange(p2) && range2->isInRange(p1)))
      return true;
    return false;
  }

  void 
  C2RPair::operator<<(const magnet::xml::Node&)
  {
    M_throw() << "Due to problems with RAll C2RPair operator<< cannot work for this class";
  }

  void 
  C2RPair::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Range") << "Pair" 
	<< magnet::xml::tag("Range1")
	<< range1
	<< magnet::xml::endtag("Range1")
	<< magnet::xml::tag("Range2")
	<< range2
	<< magnet::xml::endtag("Range2");
  }
}
