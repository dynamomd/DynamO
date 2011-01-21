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

#include "2RList.hpp"
#include "../../simulation/particle.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

C2RList::C2RList(const XMLNode& XML) 
{ operator<<(XML); }

bool 
C2RList::isInRange(const Particle&p1, const Particle&p2) const
{
  std::map<unsigned long, std::list<unsigned long> >::const_iterator
    iPtr;
  
  unsigned long a = p1.getID();
  unsigned long b = p2.getID();
  
  if (a < b)
    {
      if ((iPtr = pairmap.find(a)) != pairmap.end())
	BOOST_FOREACH(unsigned long val, iPtr->second)
	  if (b == val)
	    return true;
    }
  else
    if ((iPtr = pairmap.find(b))  != pairmap.end())
      BOOST_FOREACH(unsigned long val, iPtr->second)
	if (a == val)
	  return true;
  
  return false;
}

void 
C2RList::addPair(unsigned long a, unsigned long b)
{
  if (a < b)
    pairmap[a].push_back(b);
  else 
    pairmap[b].push_back(a);
}

const std::map<unsigned long, std::list<unsigned long> >& 
C2RList::getPairMap() const
{
  return pairmap;
}


void 
C2RList::operator<<(const XMLNode& XML)
{
  if (strcmp(XML.getAttribute("Range"),"List"))
    M_throw() << "Attempting to load a List from a non List";    
  
  try 
    {
      XMLNode xSubNode;
      
      for (long i=0; i < XML.nChildNode("RangePair"); i++)
	{ 
	  xSubNode = XML.getChildNode("RangePair",i);
	  addPair(boost::lexical_cast<unsigned long>(xSubNode.getAttribute("ID1")), 
		  boost::lexical_cast<unsigned long>(xSubNode.getAttribute("ID2")));
	}
    }
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in C2RList";
    }
}

void 
C2RList::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Range") << "List";
  typedef const std::pair<const unsigned long, std::list<unsigned long> >& thepair;
  BOOST_FOREACH(thepair mypair, pairmap)
    BOOST_FOREACH(unsigned long val, mypair.second)
    XML << xml::tag("RangePair") << xml::attr("ID1") << mypair.first
	<< xml::attr("ID2") << val << xml::endtag("RangePair");
}
