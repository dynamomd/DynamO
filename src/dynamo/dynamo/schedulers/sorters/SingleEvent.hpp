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

#include <dynamo/schedulers/sorters/event.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
  //! There is a trick used here to speed up comparisons between
  //! MinMaxHeaps.  The top element is set to HUGE_VAL, whenever the
  //! queue is cleared, or pop'd empty. This means no conditional logic
  //! is required to deal with the comparison of empty queues.
  class PELSingleEvent
  {
    Event _event;

  public:
    PELSingleEvent() { clear(); }

    inline size_t size() const { return _event.type != NONE; }
    inline bool empty() const { return _event.type == NONE; }
    inline bool full() const { return _event.type != NONE; }

    inline const Event& front() const { return _event; }
    inline const Event& top() const { return _event; }  

    inline void pop() 
    { 
      if (empty()) return;
      
      _event.type = VIRTUAL; //Force a recalculation for any further events
    }

    inline void clear() { _event.dt = HUGE_VAL; _event.type = NONE; }

    inline bool operator> (const PELSingleEvent& ip) const throw()
    { 
      return _event.dt > ip._event.dt; 
    }

    inline bool operator< (const PELSingleEvent& ip) const throw()
    { 
      return _event.dt < ip._event.dt; 
    }

    inline double getdt() const 
    {
      return _event.dt;
    }
  
    inline void stream(const double& ndt) throw()
    { _event.dt -= ndt; }

    inline void addTime(const double& ndt) throw()
    { _event.dt += ndt; }

    inline void push(const Event& __x)
    { if (__x < _event) _event = __x; }

    inline void rescaleTimes(const double& scale) throw()
    { _event.dt *= scale; }

    inline void swap(PELSingleEvent& rhs)
    { std::swap(_event, rhs._event); }  
  };
}

namespace std
{
  template<>
  inline void swap(dynamo::PELSingleEvent& lhs, dynamo::PELSingleEvent& rhs)
  { lhs.swap(rhs); }
}
