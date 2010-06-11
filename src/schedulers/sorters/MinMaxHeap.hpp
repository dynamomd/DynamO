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

#pragma once

#include "datastruct.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../datatypes/MinMaxHeap.hpp"

template<size_t Size>
class MinMaxHeapPList
{
  MinMaxHeap<intPart, Size> _innerHeap;

public:
  inline size_t size() const { return _innerHeap.size(); }
  inline bool empty() const { return _innerHeap.empty(); }
  inline bool full() const { return _innerHeap.full(); }

  inline const intPart& front() const { return _innerHeap.top(); }
  inline const intPart& top() const { return _innerHeap.top(); }  

  inline void pop() { _innerHeap.pop(); }

  inline void clear() { _innerHeap.clear(); }

  inline bool operator> (const MinMaxHeapPList& ip) const throw()
  { 
    //If the other is empty this can never be longer
    //If this is empty and the other isn't its always longer
    //Otherwise compare
    return (ip._innerHeap.empty()) 
      ? false
      : (empty() || (_innerHeap.top().dt > ip._innerHeap.top().dt)); 
  }

  inline bool operator< (const MinMaxHeapPList& ip) const throw()
  { 
    //If this is empty it can never be shorter
    //If the other is empty its always shorter
    //Otherwise compare
    return (empty()) 
      ? false 
      : (ip._innerHeap.empty() || (_innerHeap.top().dt < ip._innerHeap.top().dt)); 
  }

  inline Iflt getdt() const 
  { 
    return (empty()) ? HUGE_VAL : _innerHeap.top().dt; 
  }
  
  inline void stream(const Iflt& ndt) throw()
  {
    BOOST_FOREACH(intPart& dat, _innerHeap)
      dat.dt -= ndt;
  }

  inline void addTime(const Iflt& ndt) throw()
  {
    BOOST_FOREACH(intPart& dat, _innerHeap)
      dat.dt += ndt;
  }

  inline void push(const intPart& __x)
  {
    if (!_innerHeap.full()) { _innerHeap.insert(__x); return; }

    if (__x < _innerHeap.bottom())
      _innerHeap.replaceMax(__x);
	
    _innerHeap.unsafe_bottom().type = VIRTUAL;
  }

  inline void rescaleTimes(const Iflt& scale) throw()
  { 
    BOOST_FOREACH(intPart& dat, _innerHeap)
      dat.dt *= scale;
  }

  inline void swap(MinMaxHeapPList& rhs)
  {
    _innerHeap.swap(rhs._innerHeap);
  }
  
};

namespace std
{
  /*! \brief Template specialisation of the std::swap function for pList*/
  template<size_t Size> inline
  void swap(MinMaxHeapPList<Size>& lhs, MinMaxHeapPList<Size>& rhs)
  {
    lhs.swap(rhs);
  }
}
