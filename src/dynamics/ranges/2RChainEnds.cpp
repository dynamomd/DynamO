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

#include "2RChainEnds.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../../simulation/particle.hpp"

C2RChainEnds::C2RChainEnds(const XMLNode& XML, const DYNAMO::SimData*):
  rangeStart(0),rangeEnd(0), interval(0) 
{ 
  if (strcmp(XML.getAttribute("Range"),"ChainEnds"))
    D_throw() << "Attempting to load a ChainEnds from a "
	      << XML.getAttribute("Range");
  
  rangeStart = boost::lexical_cast<size_t>(XML.getAttribute("Start"));
  rangeEnd = boost::lexical_cast<size_t>(XML.getAttribute("End"));
  interval = boost::lexical_cast<size_t>(XML.getAttribute("Interval"));

  //Guarrantee that they are ordered
  if (rangeStart > rangeEnd)
    std::swap(rangeStart, rangeEnd);
  
  if ((rangeEnd - rangeStart + 1) % interval)
    D_throw() << "Length of range does not split into an integer"
	      << " number of intervals";
}

C2RChainEnds::C2RChainEnds(size_t r1, size_t r2, 
			   size_t l):
  rangeStart(r1),rangeEnd(r2), interval(l) 
{
  //Guarrantee that they are ordered
  if (rangeStart > rangeEnd)
    std::swap(rangeStart, rangeEnd);

  if ((rangeEnd - rangeStart + 1) % interval)
    D_throw() << "Length of range does not split into an integer"
	      << " number of intervals";
}

bool 
C2RChainEnds::isInRange(const CParticle&p1, const CParticle&p2) const
{
  if (p1.getID() > p2.getID())
    return ((p1.getID() <= rangeEnd) && (p2.getID() >= rangeStart)
	    && !((p2.getID() - rangeStart) % interval)
	    && (p1.getID() - p2.getID() == interval));
  else
    return ((p2.getID() <= rangeEnd) && (p1.getID() >= rangeStart)
	    && !((p1.getID() - rangeStart) % interval)
	    && (p2.getID() - p1.getID() == interval));
}

void 
C2RChainEnds::operator<<(const XMLNode&)
{
  D_throw() << "Due to problems with CRAll C2RChainEnds operator<<"
    " cannot work for this class";
}

void 
C2RChainEnds::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Range") << "ChainEnds" 
      << xmlw::attr("Start")
      << rangeStart
      << xmlw::attr("End")
      << rangeEnd
      << xmlw::attr("Interval")
      << interval;
}
