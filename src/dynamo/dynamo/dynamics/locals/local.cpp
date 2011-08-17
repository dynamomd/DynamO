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
#include "local.hpp"
#include "../../simulation/particle.hpp"
#include "localEvent.hpp"
#include "../ranges/1RAll.hpp"
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

Local::Local(dynamo::SimData* tmp, const char *name):
  SimBase(tmp,name),
  range(new CRAll(tmp))
{}

Local::Local(CRange* nR, dynamo::SimData* tmp, const char *name):
  SimBase(tmp, name),
  range(nR)
{}

bool 
Local::isInteraction(const Particle &p1) const
{
  return range->isInRange(p1);
}

magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream& XML, 
			    const Local& g)
{
  g.outputXML(XML);
  return XML;
}

Local* 
Local::getClass(const magnet::xml::Node& XML, dynamo::SimData* Sim)
{
  if (!strcmp(XML.getAttribute("Type"),"Wall"))
    return new LWall(XML, Sim);
  else if (!strcmp(XML.getAttribute("Type"),"RoughWall"))
    return new LRoughWall(XML, Sim);
  else if (!strcmp(XML.getAttribute("Type"),"TriangleMesh"))
    return new LTriangleMesh(XML, Sim);
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

