/*  DYNAMO:- Event driven molecular dynamics simulator 
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

#include "include.hpp"
#include <magnet/exception.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

xml::XmlStream& operator<<(xml::XmlStream& XML, 
			   const BoundaryCondition& g)
{
  g.outputXML(XML);
  return XML;
}


BoundaryCondition* 
BoundaryCondition::getClass(const magnet::xml::Node& XML, DYNAMO::SimData* tmp)
{
  if (!std::strcmp(XML.getAttribute("Boundary"),"None")
      || !std::strcmp(XML.getAttribute("Boundary"),"Null"))
    return new BCNone(tmp);
  else if (!std::strcmp(XML.getAttribute("Boundary"),"PBC"))
    return new BCPeriodic(tmp);
  else if (!std::strcmp(XML.getAttribute("Boundary"),"NoXPBC"))
    return new BCPeriodicExceptX(tmp);
  else if (!std::strcmp(XML.getAttribute("Boundary"),"OnlyXPBC"))
    return new BCPeriodicXOnly(tmp);
  else if (!std::strcmp(XML.getAttribute("Boundary"),"LE"))
    return new BCLeesEdwards(XML,tmp);
  else 
    M_throw() << XML.getAttribute("Boundary") 
	      << ", Unknown type of rectangular boundary encountered";
}
