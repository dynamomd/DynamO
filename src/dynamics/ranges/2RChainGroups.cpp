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

#include "2RChainGroups.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../../simulation/particle.hpp"

C2RChainGroups::C2RChainGroups(const XMLNode& XML, const DYNAMO::SimData*):
  range1(0),range2(0), length(0) 
{ 
  if (strcmp(XML.getAttribute("Range"),"ChainGroups"))
    D_throw() << "Attempting to load a ChainGroups from a "
	      << XML.getAttribute("Range");
  
  range1 = boost::lexical_cast<size_t>(XML.getAttribute("Start1"));
  range2 = boost::lexical_cast<size_t>(XML.getAttribute("Start2"));
  length = boost::lexical_cast<size_t>(XML.getAttribute("Length"));

  //Guarrantee that they are ordered
  if (range1 > range2)
    std::swap(range1, range2);
}

C2RChainGroups::C2RChainGroups(size_t r1, size_t r2, 
			       size_t l):
  range1(r1),range2(r2), length(l) 
{
  //Guarrantee that they are ordered
  if (range1 > range2)
    std::swap(range1, range2);
}

bool 
C2RChainGroups::isInRange(const CParticle&p1, const CParticle&p2) const
{
  if (p1.getID() > p2.getID())
    return ((((p1.getID() >= range2) && (p1.getID() < range2 + length))
	     && ((p2.getID() >= range1) && (p2.getID() < range1 + length)))
	    && (((p1.getID() - range2) % length) 
		== ((p2.getID() - range1) % length)));
  else
    return ((((p1.getID() >= range1) && (p1.getID() < range1 + length))
	     && ((p2.getID() >= range2) && (p2.getID() < range2 + length)))
	    && (((p1.getID() - range1) % length) 
		== ((p2.getID() - range2) % length)));	    
}
void 
C2RChainGroups::operator<<(const XMLNode&)
{
  D_throw() << "Due to problems with CRAll C2RChainGroups operator<<"
    " cannot work for this class";
}

void 
C2RChainGroups::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Range") << "ChainGroups" 
      << xmlw::attr("Start1")
      << range1
      << xmlw::attr("Start2")
      << range2
      << xmlw::attr("Length")
      << length;
}
