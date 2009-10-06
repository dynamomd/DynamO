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
#include "local.hpp"
#include "../../extcode/xmlParser.h"
#include "../../extcode/xmlwriter.hpp"
#include "../../simulation/particle.hpp"
#include "localEvent.hpp"
#include "../ranges/1RAll.hpp"

CLocal::CLocal(DYNAMO::SimData* tmp, const char *name):
  SimBase(tmp,name,IC_blue),
  range(new CRAll(tmp))
{}

CLocal::CLocal(CRange* nR, DYNAMO::SimData* tmp, const char *name):
  SimBase(tmp, name,IC_blue),
  range(nR)
{}

bool 
CLocal::isInteraction(const CParticle &p1) const
{
  return range->isInRange(p1);
}

xmlw::XmlStream& operator<<(xmlw::XmlStream& XML, 
			    const CLocal& g)
{
  g.outputXML(XML);
  return XML;
}

CLocal* 
CLocal::getClass(const XMLNode &XML, DYNAMO::SimData* Sim)
{
  if (!strcmp(XML.getAttribute("Type"),"Wall"))
    return new CLWall(XML, Sim);
  else if (!strcmp(XML.getAttribute("Type"),"AndersenWall"))
    return new CLAndersenWall(XML, Sim);
  else if (!strcmp(XML.getAttribute("Type"),"DoubleWall"))
    return new CLDblWall(XML, Sim);
  else if (!strcmp(XML.getAttribute("Type"),"OscillatingPlate"))
    return new CLOscillatingPlate(XML, Sim);
  else if (!strcmp(XML.getAttribute("Type"),"CylinderWall"))
    return new CLCylinder(XML, Sim);
  else 
    D_throw() << "Unknown type of Local Interaction encountered";
}

