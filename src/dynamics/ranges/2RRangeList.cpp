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

#include "2RRangeList.hpp"
#include "../../simulation/particle.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include <boost/foreach.hpp>

C2RRangeList::C2RRangeList(const XMLNode& XML, const DYNAMO::SimData* nSim):
  SimBase_const(nSim,"C2RRangeList", IC_red)
{ operator<<(XML); }

bool 
C2RRangeList::isInRange(const Particle&p1, const Particle&p2) const
{
  BOOST_FOREACH(const ClonePtr<C2Range>& rPtr, ranges)
    if (rPtr->isInRange(p1,p2))
      return true;
  
  return false;
}

void 
C2RRangeList::operator<<(const XMLNode& XML)
{
  if (strcmp(XML.getAttribute("Range"),"RangeList"))
    D_throw() << "Attempting to load a List from a non List";    
  
  try 
    {
      XMLNode xSubNode;
      
      for (long i=0; i < XML.nChildNode("RangeListItem"); i++)
	ranges.push_back(ClonePtr<C2Range>(C2Range::loadClass(XML.getChildNode("RangeListItem",i), Sim)));
    }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in C2RRangeList";
    }
}

void 
C2RRangeList::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Range") << "RangeList";

  BOOST_FOREACH(const ClonePtr<C2Range>& rPtr, ranges)
    XML << xml::tag("RangeListItem") << rPtr << xml::endtag("RangeListItem");
}
