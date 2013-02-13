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
#include <dynamo/ranges/IDRange.hpp>
#include <dynamo/ranges/IDPairRange.hpp>

namespace dynamo {
  class IDPairRangePair:public IDPairRange
  {
  public:
    IDPairRangePair(IDRange* r1, IDRange* r2 ):range1(r1), range2(r2) {}

    IDPairRangePair(const magnet::xml::Node& XML, const dynamo::Simulation* Sim)
    { 
      magnet::xml::Node subRangeXML = XML.getNode("IDRange");
      range1 = shared_ptr<IDRange>(IDRange::getClass(subRangeXML, Sim));
      ++subRangeXML;
      range2 = shared_ptr<IDRange>(IDRange::getClass(subRangeXML, Sim));
    }

    virtual bool isInRange(const Particle&p1, const Particle&p2) const
    {
      if ((range1->isInRange(p1) && range2->isInRange(p2))
	  || (range1->isInRange(p2) && range2->isInRange(p1)))
	return true;
      return false;
    }

  protected:

    virtual void outputXML(magnet::xml::XmlStream& XML) const
    {
      XML << magnet::xml::attr("Type") << "Pair"
	  << range1 << range2;
    }

    shared_ptr<IDRange> range1;
    shared_ptr<IDRange> range2;
  };
}
