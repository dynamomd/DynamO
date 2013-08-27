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
#include <dynamo/particle.hpp>

namespace dynamo {
  class IDPairRangeChains:public IDPairRange
  {
  public:
    IDPairRangeChains(unsigned long r1, unsigned long r2, unsigned long r3):
      range1(r1),range2(r2),interval(r3) 
    {
      if ((r2-r1 + 1) % r3)
	M_throw() << "Range of IDPairRangeChains does not split evenly into interval";
    }

    IDPairRangeChains(const magnet::xml::Node& XML, const dynamo::Simulation*):
      range1(0),range2(0), interval(0)
    { 
      range1 = XML.getAttribute("Start").as<unsigned long>();
      range2 = XML.getAttribute("End").as<unsigned long>();
      interval = XML.getAttribute("Interval").as<unsigned long>();
      if ((range2-range1 + 1) % interval)
	M_throw() << "Range of IDPairRangeChains does not split evenly into interval";

    }

    virtual bool isInRange(const Particle&p1, const Particle&p2) const
    {
      size_t a = std::min(p1.getID(), p2.getID()), b = std::max(p1.getID(), p2.getID());
      
      return (b - a == 1) //Check the particles are next to each other
	&& ((a >= range1) && (b <= range2)) //and within the overall range
	&& (((a  - range1) / interval) == ((b  - range1) / interval)); //and belong to the same segment
    }
  
    virtual bool isInRange(const Particle& p1) const { return (p1.getID() >= range1) && (p1.getID() <= range2); }

  protected:
    virtual void outputXML(magnet::xml::XmlStream& XML) const
    {
      XML << magnet::xml::attr("Type") << "Chains" 
	  << magnet::xml::attr("Start")
	  << range1
	  << magnet::xml::attr("End")
	  << range2
	  << magnet::xml::attr("Interval")
	  << interval;
    }

    size_t range1;
    size_t range2;
    size_t interval;
  };
}
