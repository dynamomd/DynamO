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

#include "include.hpp"
#include <magnet/exception.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <cstring>

magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream& XML, 
			   const BoundaryCondition& g)
{
  g.outputXML(XML);
  return XML;
}


BoundaryCondition* 
BoundaryCondition::getClass(const magnet::xml::Node& XML, dynamo::SimData* tmp)
{
  if (!std::strcmp(XML.getAttribute("Type"),"None")
      || !std::strcmp(XML.getAttribute("Type"),"Null"))
    return new BCNone(tmp);
  else if (!std::strcmp(XML.getAttribute("Type"),"PBC"))
    return new BCPeriodic(tmp);
  else if (!std::strcmp(XML.getAttribute("Type"),"NoXPBC"))
    return new BCPeriodicExceptX(tmp);
  else if (!std::strcmp(XML.getAttribute("Type"),"OnlyXPBC"))
    return new BCPeriodicXOnly(tmp);
  else if (!std::strcmp(XML.getAttribute("Type"),"LE"))
    return new BCLeesEdwards(XML,tmp);
  else 
    M_throw() << XML.getAttribute("Type") 
	      << ", Unknown type of boundary encountered";
}
