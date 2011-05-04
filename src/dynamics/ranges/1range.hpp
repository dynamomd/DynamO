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

#include <boost/range.hpp>
#include <iterator>

class Particle;
namespace magnet { namespace xml { class Node; } }
namespace xml { class XmlStream; }
namespace dynamo { class SimData; }
namespace { class RangeIterator; }

class CRange
{
public:
  typedef RangeIterator iterator;
  typedef RangeIterator const_iterator;

  virtual ~CRange() {};

  virtual bool isInRange(const Particle&) const = 0;

  virtual void operator<<(const magnet::xml::Node&) = 0;  

  virtual CRange* Clone() const = 0;

  virtual unsigned long size() const = 0;

  virtual unsigned long operator[](unsigned long) const = 0;

  virtual unsigned long at(unsigned long) const = 0;

  static CRange* getClass(const magnet::xml::Node&, const dynamo::SimData * Sim);

  friend xml::XmlStream& operator<<(xml::XmlStream&, const CRange&);

  virtual iterator begin() const = 0;

  virtual iterator end() const = 0;

  virtual const unsigned long& getIteratorID(const unsigned long &) const =0;// { return i; }

protected:

  virtual void outputXML(xml::XmlStream& ) const = 0;    
};

namespace {
  class RangeIterator
  {
  public:
    RangeIterator(unsigned long nPos, const CRange* nRangePtr):
      pos(nPos), rangePtr(nRangePtr) {}

    inline bool operator==(const RangeIterator& nIT) const
    { return nIT.pos == pos; }

    inline bool operator!=(const RangeIterator& nIT) const
    { return nIT.pos != pos; }

    inline RangeIterator operator+(const unsigned long& i) const
    { return RangeIterator(pos+i, rangePtr); }

    inline RangeIterator operator-(const unsigned long& i) const
    { return RangeIterator(pos-i, rangePtr); }
    
    inline RangeIterator& operator++()
    { pos++; return *this; }

    inline RangeIterator& operator++(int)
    { pos++; return *this; }

    inline const unsigned long& operator*() const
    { return rangePtr->getIteratorID(pos); }

    typedef unsigned long difference_type;
    typedef const unsigned long& value_type;
    typedef const unsigned long& reference;
    typedef const unsigned long* pointer;
    typedef std::bidirectional_iterator_tag iterator_category;

  private:
    unsigned long pos;
    const CRange* rangePtr;
  };
}
