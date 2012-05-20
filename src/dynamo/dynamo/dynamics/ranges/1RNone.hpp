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
#include <dynamo/dynamics/ranges/1range.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
  class RNone: public Range
  {
  public:
    RNone() {}

    RNone(const magnet::xml::Node& XML)
    { operator<<(XML); }

    virtual bool isInRange(const Particle&) const
    { return false; }

    virtual void operator<<(const magnet::xml::Node& XML)
    {
      if (strcmp(XML.getAttribute("Range"),"None"))
	M_throw() << "Attempting to load RNone from non None type";
    }

    virtual unsigned long size() const { return 0; }

    virtual iterator begin() const { return Range::iterator(0, this); }

    virtual iterator end() const { return Range::iterator(0, this); }

    virtual unsigned long operator[](unsigned long i) const  
    {
      M_throw() << "Nothing to access";
    }

    virtual unsigned long at(unsigned long i) const 
    { 
      M_throw() << "Nothing to access";
    }

  protected:

    virtual const unsigned long& getIteratorID(const unsigned long &i) const 
    { M_throw() << "Nothing here!"; }

    virtual void outputXML(magnet::xml::XmlStream& XML) const
    { XML << magnet::xml::attr("Range") << "All"; }
  };
}
