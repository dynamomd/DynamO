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

#include <dynamo/schedulers/sorters/event.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/containers/MinMaxHeap.hpp>

namespace dynamo {
  /*! A MinMax heap used for Particle Event Lists

    There is a trick used here to speed up comparisons between
    MinMaxHeaps.  The top element is set to HUGE_VAL, whenever the
    queue is cleared, or pop'd empty. This means no conditional logic
    is required to deal with the comparison of empty queues.
  */
  template<size_t Size>
  class PELMinMax: public magnet::containers::MinMaxHeap<Event,Size>
  {
    typedef magnet::containers::MinMaxHeap<Event,Size> Base;
  public:
    PELMinMax() { 
      clear(); 
    }

    inline void pop() { 
      Base::pop();
      if (Base::empty()) clear(); 
    }

    inline void clear() { 
      Base::clear(); 
      Base::begin()->dt = HUGE_VAL; 
    }

    inline bool operator> (const PELMinMax& ip) const {  
      return Base::begin()->dt > ip.begin()->dt;
    }

    inline bool operator< (const PELMinMax& ip) const { 
      return Base::begin()->dt < ip.begin()->dt; 
    }
    
    inline double getdt() const {
      return Base::begin()->dt;
    }
  
    inline void stream(const double& ndt) {
      for(Event& dat : *this)
	dat.dt -= ndt;
    }

    inline void push(const Event& __x) {
      if (!Base::full())
	Base::insert(__x);
      else 
	{
	  if (__x < Base::bottom())
	    Base::replaceMax(__x);
	  Base::unsafe_bottom().type = RECALCULATE;
	}
    }

    inline void rescaleTimes(const double& scale) { 
      for (Event& dat : *this)
	dat.dt *= scale;
    }

    inline void swap(PELMinMax& rhs) {
      Base::swap(rhs);
    }
  };
}

namespace std
{
  /*! \brief Template specialisation of the std::swap function for PELHeap*/
  template<size_t Size>
  void swap(dynamo::PELMinMax<Size>& lhs, dynamo::PELMinMax<Size>& rhs)
  {
    lhs.swap(rhs);
  }
}
