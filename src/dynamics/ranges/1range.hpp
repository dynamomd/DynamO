/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef CRange_H
#define CRange_H

#include <boost/range.hpp>
#include <iterator>

class XMLNode;
namespace xmlw
{
  class XmlStream;
}
class CParticle;
namespace DYNAMO
{
  class SimData;
}

class CRange
{
public:
  class iterator
  {
  public:
    iterator(unsigned long nPos, const CRange* nRangePtr):
      pos(nPos), rangePtr(nRangePtr) {}

    inline bool operator==(const iterator& nIT) const
    { return nIT.pos == pos; }

    inline bool operator!=(const iterator& nIT) const
    { return nIT.pos != pos; }

    inline iterator operator+(const unsigned long& i) const
    { return iterator(pos+i, rangePtr); }

    inline iterator operator-(const unsigned long& i) const
    { return iterator(pos-i, rangePtr); }
    
    inline iterator& operator++()
    { pos++; return *this; }

    inline iterator& operator++(int)
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

  //There are only const iterators to these sequences atm
  typedef iterator const_iterator;

  virtual ~CRange() {};

  virtual bool isInRange(const CParticle&) const = 0;

  virtual void operator<<(const XMLNode&) = 0;  

  virtual CRange* Clone() const = 0;

  virtual unsigned long size() const = 0;

  virtual unsigned long operator[](unsigned long) const = 0;

  virtual unsigned long at(unsigned long) const = 0;

  static CRange* loadClass(const XMLNode&, const DYNAMO::SimData * Sim);

  friend xmlw::XmlStream& operator<<(xmlw::XmlStream&, const CRange&);

  virtual iterator begin() const = 0;

  virtual iterator end() const = 0;

protected:
  virtual const unsigned long& getIteratorID(const unsigned long &) const =0;// { return i; }

  virtual void outputXML(xmlw::XmlStream& ) const = 0;    
};

#endif
