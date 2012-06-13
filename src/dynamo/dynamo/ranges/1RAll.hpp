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
#include <dynamo/simulation.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
  class RAll: public Range, public dynamo::SimBase_const
  {
  public:
    RAll(const dynamo::Simulation* SimDat):
      SimBase_const(SimDat,"RAll"){}

    RAll(const magnet::xml::Node& XML, const dynamo::Simulation* SimDat):
      SimBase_const(SimDat, "RAll")
    { operator<<(XML); }

    virtual bool isInRange(const Particle&) const
    { return true; }

    void operator<<(const magnet::xml::Node& XML)
    {
      if (strcmp(XML.getAttribute("Range"),"All"))
	M_throw() << "Attempting to load RAll from non All type";
    }

    virtual unsigned long size() const { return Sim->particleList.size(); }

    virtual iterator begin() const { return Range::iterator(0, this); }

    virtual iterator end() const { return Range::iterator(Sim->particleList.size(), this); }

    virtual unsigned long operator[](unsigned long i) const  
    { return i; }

    virtual unsigned long at(unsigned long i) const 
    { 
      if (i >= Sim->particleList.size())
	M_throw() << "Bad array access value in range.at()";

      return i;
    }

  protected:

    virtual const unsigned long& getIteratorID(const unsigned long &i) const { return i; }

    void outputXML(magnet::xml::XmlStream& XML) const
    { XML << magnet::xml::attr("Range") << "All"; }
  };
}
