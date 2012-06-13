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

namespace dynamo {
  class C2RNone:public C2Range
  {
  public:
    C2RNone() {}

    C2RNone(const magnet::xml::Node& XML, const dynamo::Simulation*)
    { 
      if (strcmp(XML.getAttribute("Range"),"2None"))
	M_throw() << "Attempting to load a 2None from a non 2None";
    }
    
    virtual void operator<<(const magnet::xml::Node&)
    { M_throw() << "Due to problems with RAll C2RNone operator<< cannot work for this class"; }
  
    virtual bool isInRange(const Particle&, const Particle&) const
    { return false; }
  
  protected:
    virtual void outputXML(magnet::xml::XmlStream& XML) const
    { XML << magnet::xml::attr("Range") << "2None"; }
  };
}
