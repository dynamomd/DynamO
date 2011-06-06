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
#include "1range.hpp"
#include <magnet/exception.hpp>

class CRSingle: public CRange
{
public:
  CRSingle(const magnet::xml::Node&);

  CRSingle():ID(0) {}

  virtual CRange* Clone() const { return new CRSingle(*this); };

  virtual bool isInRange(const Particle &) const;

  //The data output classes
  virtual void operator<<(const magnet::xml::Node&);
  
  virtual unsigned long size() const { return 1; };

  virtual iterator begin() const { return CRange::iterator(ID, this); }

  virtual iterator end() const { return CRange::iterator(ID+1, this); }

  virtual unsigned long operator[](unsigned long) const  
  { return ID; }

  virtual unsigned long at(unsigned long i) const 
  { 
    if (i != 0)
      M_throw() << "Bad array access value in range.at()";
    
    return ID;
  }
  
protected:
  virtual const unsigned long& getIteratorID(const unsigned long &i) const { return i; }

  virtual void outputXML(xml::XmlStream&) const;

  unsigned long ID;
};
