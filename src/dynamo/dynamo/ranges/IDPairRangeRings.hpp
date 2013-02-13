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
#include <dynamo/particle.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
  class IDPairRangeRings:public IDPairRange
  {
  public:
    IDPairRangeRings(unsigned long r1, unsigned long r2, unsigned long r3):
    range1(r1),range2(r2),interval(r3) 
    {
      if ((r2 - r1 + 1) % r3)
	M_throw() << "Range of IDPairRangeRings does not split evenly into interval";
    }

    IDPairRangeRings(const magnet::xml::Node& XML, const dynamo::Simulation*):
    range1(0),range2(0), interval(0)
    { 
      range1 = XML.getAttribute("Start").as<size_t>();
      range2 = XML.getAttribute("End").as<size_t>();
      interval = XML.getAttribute("Interval").as<size_t>();

      if ((range2-range1 + 1) % interval)
	M_throw() << "Range of IDPairRangeChains does not split evenly into interval";
    }

    virtual bool isInRange(const Particle&p1, const Particle&p2) const
    {
      if (p1.getID() > p2.getID())
	{
	  if (p1.getID() - p2.getID() != 1)
	    {
	      if ((p1.getID() - p2.getID() == interval - 1)
		  && ((p2.getID() >= range1) && (p1.getID() <= range2))
		  && (((p2.getID() - range1) / interval) == ((p1.getID() - range1) / interval)))
		return true;
	    }
	  else
	    if ((p2.getID() >= range1) && (p1.getID() <= range2))
	      if (((p2.getID() - range1) / interval) == ((p1.getID() - range1) / interval))
		return true;
	}
      else 
	{
	  if (p2.getID() - p1.getID() != 1)
	    {
	      if ((p2.getID() - p1.getID() == interval - 1)
		  && ((p1.getID() >= range1) && (p2.getID() <= range2))
		  && (((p1.getID() - range1) / interval) == ((p2.getID() - range1) / interval)))
		return true;
	    }
	  else
	    if ((p1.getID() >= range1) && (p2.getID() <= range2))
	      if (((p1.getID() - range1) / interval) == ((p2.getID() - range1) / interval))
		return true;
	}
      return false;
    }

  protected:
    virtual void outputXML(magnet::xml::XmlStream& XML) const
    {
      XML << magnet::xml::attr("Type") << "Rings" 
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
