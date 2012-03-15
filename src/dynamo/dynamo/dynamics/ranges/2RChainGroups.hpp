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
#include <dynamo/dynamics/ranges/1range.hpp>
#include <dynamo/dynamics/ranges/2range.hpp>
#include <dynamo/simulation/particle.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
  class C2RChainGroups:public C2Range
  {
  public:
    C2RChainGroups(const magnet::xml::Node& XML, const dynamo::SimData*):
      range1(0),range2(0), length(0) 
    { 
      if (strcmp(XML.getAttribute("Range"),"ChainGroups"))
	M_throw() << "Attempting to load a ChainGroups from a "
		  << XML.getAttribute("Range");
  
      range1 = XML.getAttribute("Start1").as<size_t>();
      range2 = XML.getAttribute("Start2").as<size_t>();
      length = XML.getAttribute("Length").as<size_t>();

      //Guarrantee that they are ordered
      if (range1 > range2)
	std::swap(range1, range2);
    }

    C2RChainGroups(size_t r1, size_t r2, size_t l):
      range1(r1),range2(r2), length(l) 
    {
      //Guarrantee that they are ordered
      if (range1 > range2)
	std::swap(range1, range2);
    }
  
    virtual bool isInRange(const Particle&p1, const Particle&p2) const
    {
      if (p1.getID() > p2.getID())
	return ((((p1.getID() >= range2) && (p1.getID() < range2 + length))
		 && ((p2.getID() >= range1) && (p2.getID() < range1 + length)))
		&& (((p1.getID() - range2) % length) 
		    == ((p2.getID() - range1) % length)));
      else
	return ((((p1.getID() >= range1) && (p1.getID() < range1 + length))
		 && ((p2.getID() >= range2) && (p2.getID() < range2 + length)))
		&& (((p1.getID() - range1) % length) 
		    == ((p2.getID() - range2) % length)));	    
    }
    
    virtual void operator<<(const magnet::xml::Node&)
    {
      M_throw() << "Due to problems with RAll C2RChainGroups operator<<"
	" cannot work for this class";
    }

  protected:
    virtual void outputXML(magnet::xml::XmlStream& XML) const
    {
      XML << magnet::xml::attr("Range") << "ChainGroups" 
	  << magnet::xml::attr("Start1")
	  << range1
	  << magnet::xml::attr("Start2")
	  << range2
	  << magnet::xml::attr("Length")
	  << length;
    }

    size_t range1;
    size_t range2;
    size_t length;
  };
}
