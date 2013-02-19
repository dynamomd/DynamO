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
  typedef std::vector<Event> qType;
  typedef std::priority_queue<Event, qType, 
			      std::greater<Event> > PELHeap_q_type;
  class PELHeap: public PELHeap_q_type
  {
  public:
    typedef qType::iterator CRanIt;
    typedef qType::iterator iterator;
    typedef qType::const_iterator const_iterator;
    typedef std::iterator_traits<CRanIt>::difference_type CDiff;
  
    inline iterator begin() { return c.begin(); }
    inline const_iterator begin() const { return c.begin(); }
    inline iterator end() { return c.end(); }
    inline const_iterator end() const { return c.end(); }

    inline void clear()
    { c.clear(); }

    inline bool operator> (const PELHeap& ip) const throw()
    { 
      //If the other is empty this can never be longer
      //If this is empty and the other isn't its always longer
      //Otherwise compare
      return (ip.c.empty()) 
	? false
	: (empty() || (c.front().dt > ip.c.front().dt)); 
    }

    inline bool operator< (const PELHeap& ip) const throw()
    { 
      //If this is empty it can never be shorter
      //If the other is empty its always shorter
      //Otherwise compare
      return (empty()) 
	? false 
	: (ip.c.empty() || (c.front().dt < ip.c.front().dt)); 
    }

    inline double getdt() const 
    { 
      return (c.empty()) ? HUGE_VAL : c.front().dt; 
    }
  
    inline void stream(const double& ndt) throw()
    {
      BOOST_FOREACH(Event& dat, c)
	dat.dt -= ndt;
    }

    inline void push(const Event& __x)
    {
      c.push_back(__x);
      std::push_heap(c.begin(), c.end(), comp);
    }

    inline void rescaleTimes(const double& scale) throw()
    { 
      BOOST_FOREACH(Event& dat, c)
	dat.dt *= scale;
    }

    inline void swap(PELHeap& rhs)
    {
      c.swap(rhs.c);
    }
  };
}

namespace std
{
  /*! \brief Template specialisation of the std::swap function for pList*/
  template<> inline void swap(dynamo::PELHeap& lhs, dynamo::PELHeap& rhs)
  { lhs.swap(rhs); }
}
