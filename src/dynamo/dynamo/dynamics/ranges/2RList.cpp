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

#include "2RList.hpp"
#include "../../simulation/particle.hpp"
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>

C2RList::C2RList(const magnet::xml::Node& XML) 
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
	{
	  BOOST_FOREACH(unsigned long val, iPtr->second)
	    if (b == val)
	      return true;
	}
    }
  else
    if ((iPtr = pairmap.find(b))  != pairmap.end())
      {
	BOOST_FOREACH(unsigned long val, iPtr->second)
	  if (a == val)
	    return true;
      }
  
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
C2RList::operator<<(const magnet::xml::Node& XML)
{
  if (strcmp(XML.getAttribute("Range"),"List"))
    M_throw() << "Attempting to load a List from a non List";    
  
  try 
    {
      for (magnet::xml::Node node = XML.getNode("RangePair"); node.valid(); ++node)
	addPair(node.getAttribute("ID1").as<unsigned long>(), 
		node.getAttribute("ID2").as<unsigned long>());
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
