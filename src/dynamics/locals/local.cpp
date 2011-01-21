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
#include "local.hpp"
#include "../../extcode/xmlParser.h"
#include "../../extcode/xmlwriter.hpp"
#include "../../simulation/particle.hpp"
#include "localEvent.hpp"
#include "../ranges/1RAll.hpp"

Local::Local(DYNAMO::SimData* tmp, const char *name):
  SimBase(tmp,name,IC_blue),
  range(new CRAll(tmp))
{}

Local::Local(CRange* nR, DYNAMO::SimData* tmp, const char *name):
  SimBase(tmp, name,IC_blue),
  range(nR)
{}

bool 
Local::isInteraction(const Particle &p1) const
{
  return range->isInRange(p1);
}

xml::XmlStream& operator<<(xml::XmlStream& XML, 
			    const Local& g)
{
  g.outputXML(XML);
  return XML;
}

Local* 
Local::getClass(const XMLNode &XML, DYNAMO::SimData* Sim)
{
  if (!strcmp(XML.getAttribute("Type"),"Wall"))
    return new CLWall(XML, Sim);
  else if (!strcmp(XML.getAttribute("Type"),"RoughWall"))
    return new LRoughWall(XML, Sim);
  else if (!strcmp(XML.getAttribute("Type"),"AndersenWall"))
    return new CLAndersenWall(XML, Sim);
  else if (!strcmp(XML.getAttribute("Type"),"DoubleWall"))
    return new CLDblWall(XML, Sim);
  else if (!strcmp(XML.getAttribute("Type"),"OscillatingPlate"))
    return new CLOscillatingPlate(XML, Sim);
  else if (!strcmp(XML.getAttribute("Type"),"CylinderWall"))
    return new CLCylinder(XML, Sim);
  else if (!strcmp(XML.getAttribute("Type"),"SphereWall"))
    return new CLSphere(XML, Sim);
  else 
    M_throw() << XML.getAttribute("Type")
	      << ", Unknown type of Local Interaction encountered";
}

