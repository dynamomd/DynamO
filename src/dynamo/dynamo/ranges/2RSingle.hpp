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
#include <dynamo/ranges/1range.hpp>
#include <dynamo/ranges/2range.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
  class C2RSingle:public C2Range
  {
  public:
    C2RSingle(Range* r1):range(r1) {}
  
    C2RSingle(const magnet::xml::Node& XML, const dynamo::SimData* Sim)
    { 
      if (strcmp(XML.getAttribute("Range"),"2Single"))
	M_throw() << "Attempting to load a 2Single from a non pair";
  
      range = shared_ptr<Range>(Range::getClass(XML.getNode("SingleRange"), Sim));
    }

    virtual bool isInRange(const Particle&p1, const Particle&p2) const
    {
      return (range->isInRange(p1) && range->isInRange(p2));
    }

    virtual void operator<<(const magnet::xml::Node&)
    {
      M_throw() << "Due to problems with C2RSingle operator<< cannot work for this class";
    }
  
    const shared_ptr<Range>& getRange() const { return range; }

  protected:
    virtual void outputXML(magnet::xml::XmlStream& XML) const
    {
      XML << magnet::xml::attr("Range") << "2Single" 
	  << magnet::xml::tag("SingleRange")
	  << range
	  << magnet::xml::endtag("SingleRange");
    }

    shared_ptr<Range> range;
  };
}
