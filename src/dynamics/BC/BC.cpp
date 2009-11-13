/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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
#include "../../extcode/xmlParser.h"
#include "../../extcode/xmlwriter.hpp"
#include "../../base/is_exception.hpp"
#include <string.h>

xmlw::XmlStream& operator<<(xmlw::XmlStream& XML, 
			    const BoundaryCondition& g)
{
  g.outputXML(XML);
  return XML;
}


BoundaryCondition* 
BoundaryCondition::loadClass(const XMLNode &XML, DYNAMO::SimData* tmp)
{
  if (!strcmp(XML.getAttribute("Boundary"),"None")
      || !strcmp(XML.getAttribute("Boundary"),"Null"))
    return new BCNone(tmp);
  else if (!strcmp(XML.getAttribute("Shape"),"Square"))
    {
      if (!strcmp(XML.getAttribute("Boundary"),"PBC"))
	return new BCSquarePeriodic(tmp);
      else if (!strcmp(XML.getAttribute("Boundary"),"LE"))
	return new BCSquareLeesEdwards(XML,tmp);
      else 
	D_throw() << "Unknown type of square boundary encountered";
    } 
  else if (!strcmp(XML.getAttribute("Shape"),"Rectangular"))
    {
      if (!strcmp(XML.getAttribute("Boundary"),"PBC"))
	return new BCRectangularPeriodic(tmp);
      else if (!strcmp(XML.getAttribute("Boundary"),"NoXPBC"))
	return new BCSquarePeriodicExceptX(tmp);
      else if (!strcmp(XML.getAttribute("Boundary"),"OnlyXPBC"))
	return new BCSquarePeriodicXOnly(tmp);
      else if (!strcmp(XML.getAttribute("Boundary"),"LE"))
	return new BCRectangularLeesEdwards(XML,tmp);
      else 
	D_throw() << "Unknown type of rectangular boundary encountered";
    }
  else
    D_throw() << "Unknown shape of Boundary encountered";

}
