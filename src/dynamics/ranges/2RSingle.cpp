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

#include "2RSingle.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"

C2RSingle::C2RSingle(const XMLNode& XML, const DYNAMO::SimData* Sim):
  range(NULL)
{ 

  if (strcmp(XML.getAttribute("Range"),"2Single"))
    I_throw() << "Attempting to load a 2Single from a non pair";
  
  range.set_ptr(CRange::loadClass(XML.getChildNode("SingleRange"), Sim));
}

bool 
C2RSingle::isInRange(const CParticle&p1, const CParticle&p2) const
{
  return (range->isInRange(p1) && range->isInRange(p2));
}

void 
C2RSingle::operator<<(const XMLNode&)
{
  I_throw() << "Due to problems with C2RSingle operator<< cannot work for this class";
}

void 
C2RSingle::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Range") << "2Single" 
      << xmlw::tag("SingleRange")
      << range
      << xmlw::endtag("SingleRange");
}

