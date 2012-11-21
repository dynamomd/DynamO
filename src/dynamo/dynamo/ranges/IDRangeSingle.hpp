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
#include <dynamo/particle.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
  class IDRangeSingle: public IDRange
  {
  public:
    IDRangeSingle(const magnet::xml::Node& XML) 
    { operator<<(XML); }

    IDRangeSingle():ID(0) {}

    virtual bool isInRange(const Particle &part) const
    { return part.getID() == ID; }

    virtual void operator<<(const magnet::xml::Node& XML)
    {
      if (strcmp(XML.getAttribute("Range"),"Single"))
	M_throw() << "Attempting to load IDRangeSingle from non single";
      try {
	ID = XML.getAttribute("ID").as<size_t>();
      }
      catch (boost::bad_lexical_cast &)
	{ M_throw() << "Failed a lexical cast in IDRangeSingle"; }
    }

    virtual unsigned long size() const { return 1; };

    virtual unsigned long operator[](unsigned long) const  
    { return ID; }

    virtual unsigned long at(unsigned long i) const 
    { 
      if (i != 0)
	M_throw() << "Bad array access value in range.at()";
    
      return ID;
    }
  
  protected:

    virtual void outputXML(magnet::xml::XmlStream& XML) const
    {
      XML << magnet::xml::attr("Range") << "Single"
	  << magnet::xml::attr("ID") << ID;
    }

    unsigned long ID;
  };
}
