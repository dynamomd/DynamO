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

#pragma once
#include <dynamo/dynamics/ranges/2range.hpp>
#include <dynamo/simulation/particle.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>
#include <boost/unordered_set.hpp>

namespace dynamo {
  class C2RList:public C2Range
  {
    typedef std::pair<unsigned long, unsigned long> Key;
    typedef boost::unordered_set<Key> Container;

  public:
    C2RList(const magnet::xml::Node& XML) 
    { operator<<(XML); }

    C2RList() {}

    virtual bool isInRange(const Particle&p1, const Particle&p2) const
    {
      return pairmap.find(Key(std::min(p1.getID(), p2.getID()), std::max(p1.getID(), p2.getID()))) == pairmap.end();
    }

    void addPair(unsigned long a, unsigned long b)
    { pairmap.insert(Key(std::min(a,b), std::max(a,b))); }

    const Container& getPairMap() const { return pairmap; }

    virtual void operator<<(const magnet::xml::Node& XML)
    {
      if (strcmp(XML.getAttribute("Range"),"List"))
	M_throw() << "Attempting to load a List from a non List";    
  
      try 
	{
	  for (magnet::xml::Node node = XML.fastGetNode("RangePair"); node.valid(); ++node)
	    addPair(node.getAttribute("ID1").as<unsigned long>(), 
		    node.getAttribute("ID2").as<unsigned long>());
	}
      catch (boost::bad_lexical_cast &)
	{
	  M_throw() << "Failed a lexical cast in C2RList";
	}
    }

  protected:
    virtual void outputXML(magnet::xml::XmlStream& XML) const
    {
      XML << magnet::xml::attr("Range") << "List";
      BOOST_FOREACH(const Key& key, pairmap)
	XML << magnet::xml::tag("RangePair") 
	    << magnet::xml::attr("ID1") << key.first
	    << magnet::xml::attr("ID2") << key.second
	    << magnet::xml::endtag("RangePair");
    }

    Container pairmap;
  };
}
