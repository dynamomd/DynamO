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

#include <dynamo/dynamics/ranges/2RRing.hpp>
#include <dynamo/simulation/particle.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <cstring>

namespace dynamo {
  C2RRing::C2RRing(const magnet::xml::Node& XML, const dynamo::SimData*):
    range1(0),range2(0) 
  { 
    if (strcmp(XML.getAttribute("Range"),"Ring"))
      M_throw() << "Attempting to load a ring from a non ring";
  
    range1 = XML.getAttribute("Start").as<unsigned long>();
    range2 = XML.getAttribute("End").as<unsigned long>();
  }

  bool 
  C2RRing::isInRange(const Particle&p1, const Particle&p2) const
  {
    if (p1.getID() > p2.getID())
      {
	if ((p1.getID() - p2.getID()) == 1)
	  {
	    if ((p2.getID() >= range1) && (p1.getID() <= range2))
	      return true;
	    else
	      return false;
	  }
	else
	  if ((p1.getID() == range2) && (p2.getID() == range1))
	    return true;
	  else 
	    return false;
      }
    else
      {
	if ((p2.getID() - p1.getID()) == 1)
	  {
	    if ((p1.getID() >= range1) && (p2.getID() <= range2))
	      return true;
	    else 
	      return false;
	  }
	else
	  if ((p2.getID() == range2) && (p1.getID() == range1))
	    return true;
	  else
	    return false;
      }

    return false;
  }

  void 
  C2RRing::operator<<(const magnet::xml::Node&)
  {
    M_throw() << "Due to problems with CRAll C2RRing::operator<< cannot work for this class";
  }

  void 
  C2RRing::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Range") << "Ring" 
	<< magnet::xml::attr("Start")
	<< range1
	<< magnet::xml::attr("End")
	<< range2;
  }
}
