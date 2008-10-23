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

#ifndef CSSDataStruct_H
#define CSSDataStruct_H

#include "../../dynamics/eventtypes.hpp"
#include "../../dynamics/interactions/intEvent.hpp"
#include "../../dynamics/globals/globEvent.hpp"
#include <boost/foreach.hpp>
#include <boost/pool/pool.hpp>
#include <boost/pool/pool_alloc.hpp>
#include <queue>
#include <deque>

//Datatype for a single event, stored in lists for each particle
struct intPart
{
  inline intPart(Iflt ndt, int direction) throw():
    dt(ndt),
    p2(0),
    type(CELL),
    collCounter2(direction)
  {}
  
  inline intPart(Iflt ndt, EEventType nT, const CParticle& nID2, unsigned long long nCC2) throw():
    dt(ndt),
    p2(nID2.getID()),
    type(nT),
    collCounter2(nCC2)
  {}

  inline intPart(Iflt ndt, EEventType nT) throw():
    dt(ndt),
    p2(0),
    type(nT)
  {}

  inline intPart(const CIntEvent& coll, unsigned long long nCC2) throw():
    dt(coll.getdt()),
    p2(coll.getParticle2().getID()),
    type(INTERACTION),
    collCounter2(nCC2)
  {    
    if (coll.getType() == NONE) type = NONE;
  }

  inline intPart(const CGlobEvent& coll) throw():
    dt(coll.getdt()),
    p2(0),
    type(GLOBAL)
  {
    if (coll.getType() == NONE) type = NONE;
  }

  inline bool operator< (const intPart& ip) const throw()
  { return dt < ip.dt; }

  inline bool operator> (const intPart& ip) const throw()
  { return dt > ip.dt; }

  inline void stream(Iflt ndt) throw() { dt -= ndt; }

  mutable Iflt dt;
  size_t p2;
  EEventType type;
  unsigned long long collCounter2;
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
  { return c.front().dt > ip.c.front().dt; }

  inline bool operator< (const pList& ip) const throw()
  { return c.front().dt < ip.c.front().dt; }

  inline const Iflt& getdt() const { return c.front().dt; }
  
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
};
#endif
