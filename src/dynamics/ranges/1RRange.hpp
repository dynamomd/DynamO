/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef C1RRange_H
#define C1RRange_H

#include "1range.hpp"
#include "../../base/is_exception.hpp"
#include "../../simulation/particle.hpp"

class CRRange: public CRange
{
public:
  CRRange(const XMLNode&);
  CRRange(unsigned int s, unsigned int e):startID(s),endID(e) {}

  virtual CRange* Clone() const { return new CRRange(*this); }

  inline virtual bool isInRange(const Particle& part) const
  {
    if ((part.getID() >= startID) && (part.getID() <= endID))
      return true;
    return false;
  }

  //The data output classes
  virtual void operator<<(const XMLNode&);
  
  virtual unsigned long size() const { return (endID - startID + 1); };

  unsigned long getStart() { return startID; }
  unsigned long getEnd() { return endID; }
  
  virtual iterator begin() const { return CRange::iterator(startID, this); }

  virtual iterator end() const { return CRange::iterator(endID+1, this); }

  virtual unsigned long operator[](unsigned long i) const  
  { return startID + i; }

  virtual unsigned long at(unsigned long i) const 
  { 
    if (i > endID - startID)
      D_throw() << "Bad array access value in range.at()";
    
    return startID + i;
  }

protected:
  virtual void outputXML(xmlw::XmlStream&) const;

  virtual const unsigned long& getIteratorID(const unsigned long &i) const { return i; }

  unsigned long startID;
  unsigned long endID;
};

#endif
