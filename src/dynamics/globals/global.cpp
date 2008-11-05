/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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
#include "global.hpp"
#include "../../extcode/xmlParser.h"
#include "../../extcode/xmlwriter.hpp"
#include "../../simulation/particle.hpp"
#include "globEvent.hpp"
#include "../ranges/1RAll.hpp"


CGlobal::CGlobal(const DYNAMO::SimData* tmp, const char *name):
  SimBase_const(tmp,name, IC_blue),
  range(new CRAll(tmp))
{}

CGlobal::CGlobal(CRange* nR, const DYNAMO::SimData* tmp, const char *name):
  SimBase_const(tmp, name, IC_blue),
  range(nR)
{}

bool 
CGlobal::isInteraction(const CParticle &p1) const
{
  return range->isInRange(p1);
}

xmlw::XmlStream& operator<<(xmlw::XmlStream& XML, 
			    const CGlobal& g)
{
  g.outputXML(XML);
  return XML;
}

CGlobal* 
CGlobal::getClass(const XMLNode &XML, const DYNAMO::SimData* Sim)
{
  if (!strcmp(XML.getAttribute("Type"),"AndersenWall"))
    return new CGAndersenWall(XML, Sim);
  else if (!strcmp(XML.getAttribute("Type"),"Wall"))
    return new CGWall(XML, Sim);
  else if (!strcmp(XML.getAttribute("Type"),"Cells"))
    return new CGCells(XML, Sim);
  else if (!strcmp(XML.getAttribute("Type"),"ListAndCell"))
    return new CGListAndCell(XML, Sim);
  else 
    D_throw() << "Unknown type of Global Interaction encountered";
}

