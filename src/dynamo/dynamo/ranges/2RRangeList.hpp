/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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
#include <dynamo/ranges/2range.hpp>
#include <dynamo/base.hpp>
#include <dynamo/particle.hpp>
#include <boost/foreach.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <list>

namespace dynamo {
  class C2RRangeList:public C2Range, dynamo::SimBase_const
  {
  public:
    C2RRangeList(const dynamo::SimData* nSim):SimBase_const(nSim,"C2RRangeList") {}

    C2RRangeList(const magnet::xml::Node& XML, const dynamo::SimData* nSim):
      SimBase_const(nSim, "C2RRangeList")
    { operator<<(XML); }
    
    virtual bool isInRange(const Particle&p1, const Particle&p2) const
    {
      BOOST_FOREACH(const shared_ptr<C2Range>& rPtr, ranges)
	if (rPtr->isInRange(p1,p2))
	  return true;
  
      return false;
    }

    virtual void operator<<(const magnet::xml::Node& XML)
    {
      if (strcmp(XML.getAttribute("Range"),"RangeList"))
	M_throw() << "Attempting to load a List from a non List";    
  
      try 
	{
	  for (magnet::xml::Node node = XML.fastGetNode("RangeListItem"); node.valid(); ++node)
	    ranges.push_back(shared_ptr<C2Range>(C2Range::getClass(node, Sim)));
	}
      catch (boost::bad_lexical_cast &)
	{
	  M_throw() << "Failed a lexical cast in C2RRangeList";
	}
    }

    void addRange(C2Range* nRange)
    { ranges.push_back(shared_ptr<C2Range>(nRange)); }
  
  protected:
    virtual void outputXML(magnet::xml::XmlStream& XML) const
    {
      XML << magnet::xml::attr("Range") << "RangeList";

      BOOST_FOREACH(const shared_ptr<C2Range>& rPtr, ranges)
	XML << magnet::xml::tag("RangeListItem") << rPtr << magnet::xml::endtag("RangeListItem");
    }

    std::list<shared_ptr<C2Range> > ranges;
  };
}
