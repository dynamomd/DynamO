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
#include <dynamo/ranges/IDPairRange.hpp>
#include <dynamo/base.hpp>
#include <dynamo/particle.hpp>
#include <boost/foreach.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <list>

namespace dynamo {
  class IDPairRangeUnion:public IDPairRange, dynamo::SimBase_const
  {
  public:
    IDPairRangeUnion(const dynamo::Simulation* nSim):
      SimBase_const(nSim, "IDPairRangeUnion") {}

    IDPairRangeUnion(const magnet::xml::Node& XML, const dynamo::Simulation* nSim):
      SimBase_const(nSim, "IDPairRangeUnion")
    {
      for (magnet::xml::Node node = XML.fastGetNode("IDPairRange"); node.valid(); ++node)
	{
	  shared_ptr<IDPairRange> ptr(IDPairRange::getClass(node, Sim));
	  ranges.push_back(ptr);
	}
    }
    
    virtual bool isInRange(const Particle&p1, const Particle&p2) const
    {
      BOOST_FOREACH(const shared_ptr<IDPairRange>& rPtr, ranges)
	if (rPtr->isInRange(p1,p2))
	  return true;
  
      return false;
    }

    void addRange(IDPairRange* nRange)
    { ranges.push_back(shared_ptr<IDPairRange>(nRange)); }
  
  protected:
    virtual void outputXML(magnet::xml::XmlStream& XML) const
    {
      XML << magnet::xml::attr("Type") << "Union";

      BOOST_FOREACH(const shared_ptr<IDPairRange>& rPtr, ranges)
	XML << rPtr;
    }

    std::list<shared_ptr<IDPairRange> > ranges;
  };
}
