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
#include <dynamo/base/is_simdata.hpp>

namespace dynamo {
  class RNone: public Range
  {
  public:
    RNone() {}

    RNone(const magnet::xml::Node&);

    virtual bool isInRange(const Particle&) const
    { return false; }

    //The data output classes
    virtual void operator<<(const magnet::xml::Node&);

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

    virtual void outputXML(magnet::xml::XmlStream&) const;
  };
}
