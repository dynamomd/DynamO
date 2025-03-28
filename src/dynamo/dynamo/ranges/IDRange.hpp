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
#include <iterator>
#include <memory>

namespace magnet {
namespace xml {
class Node;
class XmlStream;
} // namespace xml
} // namespace magnet

namespace dynamo {
using std::shared_ptr;
class Simulation;
class Particle;

class IDRange {
public:
  class iterator {
    friend class IDRange;

    iterator(unsigned long nPos, const IDRange *nRangePtr)
        : pos(nPos), rangePtr(nRangePtr) {}

  public:
    inline bool operator==(const iterator &nIT) const { return nIT.pos == pos; }

    inline bool operator!=(const iterator &nIT) const { return nIT.pos != pos; }

    inline iterator operator+(const unsigned long &i) const {
      return iterator(pos + i, rangePtr);
    }

    inline iterator operator-(const unsigned long &i) const {
      return iterator(pos - i, rangePtr);
    }

    inline iterator &operator++() {
      ++pos;
      return *this;
    }

    inline iterator &operator++(int) {
      ++pos;
      return *this;
    }

    inline size_t operator*() const { return (*rangePtr)[pos]; }

    typedef size_t difference_type;
    typedef size_t value_type;
    typedef size_t reference;
    typedef size_t *pointer;
    typedef std::bidirectional_iterator_tag iterator_category;

  private:
    size_t pos;
    const IDRange *rangePtr;
  };

  typedef iterator const_iterator;

  virtual ~IDRange() {};

  virtual bool isInRange(const Particle &) const = 0;

  virtual unsigned long size() const = 0;

  virtual unsigned long operator[](unsigned long) const = 0;

  virtual unsigned long at(unsigned long) const = 0;

  static IDRange *getClass(const magnet::xml::Node &,
                           const dynamo::Simulation *Sim);

  friend magnet::xml::XmlStream &operator<<(magnet::xml::XmlStream &XML,
                                            const IDRange &range);

  iterator begin() const { return IDRange::iterator(0, this); }
  iterator end() const { return IDRange::iterator(size(), this); }

  inline bool empty() const { return begin() == end(); }

protected:
  virtual void outputXML(magnet::xml::XmlStream &) const = 0;
};
} // namespace dynamo
