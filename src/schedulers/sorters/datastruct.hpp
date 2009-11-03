/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef CSSDataStruct_H
#define CSSDataStruct_H

#include "../../dynamics/eventtypes.hpp"
#include "../../dynamics/interactions/intEvent.hpp"
#include "../../dynamics/globals/globEvent.hpp"
#include "../../dynamics/locals/localEvent.hpp"
#include <boost/foreach.hpp>
#include <boost/pool/pool.hpp>
#include <boost/pool/pool_alloc.hpp>
#include <queue>
#include <deque>
#include "../../dynamics/globals/global.hpp"

//Datatype for a single event, stored in lists for each particle
class intPart
{
public:   
  inline intPart():
    dt(HUGE_VAL),
    collCounter2(-1),
    type(NONE),
    p2(-1)    
  {}

  inline intPart(const Iflt& ndt, const unsigned long long& direction) throw():
    dt(ndt),
    collCounter2(direction),
    type(CELL),
    p2(0)
  {}
  
  inline intPart(const Iflt& ndt, const EEventType& nT, 
		 const size_t& nID2, const unsigned long long& nCC2) throw():
    dt(ndt),
    collCounter2(nCC2),
    type(nT),
    p2(nID2)
    {}

  inline intPart(const Iflt& ndt, const EEventType& nT) throw():
    dt(ndt),
    type(nT),
    p2(0)
  {}

  inline intPart(const CIntEvent& coll, const unsigned long long& nCC2) throw():
    dt(coll.getdt()),
    collCounter2(nCC2),
    type(INTERACTION),
    p2(coll.getParticle2ID())
  {}

  inline intPart(const CGlobEvent& coll) throw():
    dt(coll.getdt()),
    type(GLOBAL),
    p2(coll.getGlobalID())
  {}

  inline intPart(const CLocalEvent& coll) throw():
    dt(coll.getdt()),
    type(LOCAL),
    p2(coll.getLocalID())
  {}

  inline bool operator< (const intPart& ip) const throw()
  { return dt < ip.dt; }

  inline bool operator> (const intPart& ip) const throw()
  { return dt > ip.dt; }

  inline void stream(const Iflt& ndt) throw() { dt -= ndt; }

  mutable Iflt dt;
  unsigned long long collCounter2;
  EEventType type;
  size_t p2;  
};

typedef std::vector<intPart> qType;
typedef std::priority_queue<intPart, qType, 
			    std::greater<intPart> > pList_q_type;
class pList: public pList_q_type
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

  inline bool operator> (const pList& ip) const throw()
  { 
    //If the other is empty this can never be longer
    //If this is empty and the other isn't its always longer
    //Otherwise compare
    return (ip.c.empty()) 
      ? false
      : (empty() || (c.front().dt > ip.c.front().dt)); 
  }

  inline bool operator< (const pList& ip) const throw()
  { 
    //If this is empty it can never be shorter
    //If the other is empty its always shorter
    //Otherwise compare
    return (empty()) 
      ? false 
      : (ip.c.empty() || (c.front().dt < ip.c.front().dt)); 
  }

  inline Iflt getdt() const 
  { 
    return (c.empty()) ? HUGE_VAL : c.front().dt; 
  }
  
  inline void stream(const Iflt& ndt) throw()
  {
    BOOST_FOREACH(intPart& dat, c)
      dat.dt -= ndt;
  }

  inline void addTime(const Iflt& ndt) throw()
  {
    BOOST_FOREACH(intPart& dat, c)
      dat.dt += ndt;
  }

  inline void push(const intPart& __x)
  {
    c.push_back(__x);
    std::push_heap(c.begin(), c.end(), comp);
  }

  inline void rescaleTimes(const Iflt& scale) throw()
  { 
    BOOST_FOREACH(intPart& dat, c)
      dat.dt *= scale;
  }

  inline void swap(pList& rhs)
  {
    c.swap(rhs.c);
  }
};

namespace std
{
  /*! \brief Template specialisation of the std::swap function for pList*/
  template<> inline
  void swap(pList& lhs, pList& rhs)
  {
    lhs.swap(rhs);
  }
}
#endif
