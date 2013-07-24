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
#include <queue>

namespace dynamo {
  class PELHeap: public std::priority_queue<Event, std::vector<Event>, std::greater<Event> >
  {
    typedef std::priority_queue<Event, std::vector<Event>, std::greater<Event> > Base;
  public:
    inline bool operator> (const PELHeap& ip) const { 
      //If the other is empty this can never be longer
      //If this is empty and the other isn't its always longer
      //Otherwise compare
      return (ip.empty()) ?  false : (Base::empty() || (Base::top().dt > ip.top().dt)); 
    }

    inline void clear() {
      c.clear();
    }
    
    inline bool operator< (const PELHeap& ip) const {
      return (ip > *this);
    }

    inline double getdt() const {
      return (c.empty()) ? HUGE_VAL : Base::top().dt; 
    }
  
    inline void stream(const double& ndt) {
      for (Event& dat : c)
	dat.dt -= ndt;
    }

    inline void rescaleTimes(const double& scale) { 
      for (Event& dat : c)
	dat.dt *= scale;
    }

    inline void swap(PELHeap& rhs) {
      std::swap(c, rhs.c);
    }
  };
}

namespace std
{
  /*! \brief Template specialisation of the std::swap function for pList*/
  template<> inline void swap(dynamo::PELHeap& lhs, dynamo::PELHeap& rhs)
  { lhs.swap(rhs); }
}
