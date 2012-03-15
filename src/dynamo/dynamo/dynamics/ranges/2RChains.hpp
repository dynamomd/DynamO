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

namespace dynamo {
  class C2RChains:public C2Range
  {
  public:
    C2RChains(unsigned long r1, unsigned long r2, unsigned long r3):
      range1(r1),range2(r2),interval(r3) 
    {
      if ((r2-r1 + 1) % r3)
	M_throw() << "Range of C2RChains does not split evenly into interval";
    }

    C2RChains(const magnet::xml::Node& XML, const dynamo::SimData*):
      range1(0),range2(0), interval(0)
    { 
      if (strcmp(XML.getAttribute("Range"),"Chains"))
	M_throw() << "Attempting to load a chains from a non chains";
  
      range1 = XML.getAttribute("Start").as<unsigned long>();
      range2 = XML.getAttribute("End").as<unsigned long>();
      interval = XML.getAttribute("Interval").as<unsigned long>();
      if ((range2-range1 + 1) % interval)
	M_throw() << "Range of C2RChains does not split evenly into interval";

    }

    virtual bool isInRange(const Particle&p1, const Particle&p2) const
    {
      size_t a = std::min(p1.getID(), p2.getID()), b = std::max(p1.getID(), p2.getID());
      return (b - a == 1)
	&& ((a >= range1) && (b <= range2))
	&& (((a  - range1) / interval) == ((b  - range1) / interval));
    }

    virtual void operator<<(const magnet::xml::Node&)
    {
      M_throw() << "Due to problems with RAll C2RChains::operator<< cannot work for this class";
    }
  
  protected:
    virtual void outputXML(magnet::xml::XmlStream& XML) const
    {
      XML << magnet::xml::attr("Range") << "Chains" 
	  << magnet::xml::attr("Start")
	  << range1
	  << magnet::xml::attr("End")
	  << range2
	  << magnet::xml::attr("Interval")
	  << interval;
    }

    unsigned long range1;
    unsigned long range2;
    unsigned long interval;
  };
}
