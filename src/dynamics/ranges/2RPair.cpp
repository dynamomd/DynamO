/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "2RPair.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"

C2RPair::C2RPair(const XMLNode& XML, const DYNAMO::SimData* Sim):
  range1(NULL),range2(NULL) 
{ 

  if (strcmp(XML.getAttribute("Range"),"Pair"))
    M_throw() << "Attempting to load a pair from a non pair";
  
  XMLNode xSubNode = XML.getChildNode("Range1");
  range1.set_ptr(CRange::loadClass(xSubNode, Sim));
  
  xSubNode = XML.getChildNode("Range2");
  range2.set_ptr(CRange::loadClass(xSubNode, Sim));
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
C2RPair::operator<<(const XMLNode&)
{
  M_throw() << "Due to problems with CRAll C2RPair operator<< cannot work for this class";
}

void 
C2RPair::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Range") << "Pair" 
      << xml::tag("Range1")
      << range1
      << xml::endtag("Range1")
      << xml::tag("Range2")
      << range2
      << xml::endtag("Range2");
}

