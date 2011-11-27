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
#include <dynamo/simulation/particle.hpp>
#include <magnet/exception.hpp>

namespace dynamo {
  class RRange: public Range
  {
  public:
    RRange(const magnet::xml::Node&);
    RRange(unsigned int s, unsigned int e):startID(s),endID(e) {}

    inline virtual bool isInRange(const Particle& part) const
    {
      if ((part.getID() >= startID) && (part.getID() <= endID))
	return true;
      return false;
    }

    //The data output classes
    virtual void operator<<(const magnet::xml::Node&);
  
    virtual unsigned long size() const { return (endID - startID + 1); };

    unsigned long getStart() { return startID; }
    unsigned long getEnd() { return endID; }
  
    virtual iterator begin() const { return Range::iterator(startID, this); }

    virtual iterator end() const { return Range::iterator(endID+1, this); }

    virtual unsigned long operator[](unsigned long i) const  
    { return startID + i; }

    virtual unsigned long at(unsigned long i) const 
    { 
      if (i > endID - startID)
	M_throw() << "Bad array access value in range.at()";
    
      return startID + i;
    }

  protected:
    virtual void outputXML(magnet::xml::XmlStream&) const;

    virtual const unsigned long& getIteratorID(const unsigned long &i) const { return i; }

    unsigned long startID;
    unsigned long endID;
  };
}
