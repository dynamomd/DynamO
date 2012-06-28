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
#include <dynamo/eventtypes.hpp>
#include <dynamo/interactions/intEvent.hpp>
#include <dynamo/globals/globEvent.hpp>
#include <dynamo/locals/localEvent.hpp>
#include <dynamo/globals/global.hpp>
#include <boost/foreach.hpp>
#include <algorithm>

namespace dynamo {
  /*! \brief A generic event type, which the more specialised events
      are converted to before they are sorted.

      This conversion is lossy, so events need to be recalculated if
      they are to be exectuted. 

      The RECALCULATE event type is special. If any IntEvent, GlobalEvent
      or LocalEvent has a type RECALCULATE, it is carried through. RECALCULATE
      events cause the system to be moved forward in time and the
      events for the particle are recalculated. This can all be
      handled by the scheduler.
   */
  class Event
  {
  public:   
    inline Event():
      dt(HUGE_VAL),
      collCounter2(std::numeric_limits<unsigned long>::max()),
      type(NONE),
      p2(std::numeric_limits<size_t>::max())    
    {}

    inline Event(const double& ndt, const EEventType& nT, 
		 const size_t& nID2, const unsigned long & nCC2) throw():
      dt(ndt),
      collCounter2(nCC2),
      type(nT),
      p2(nID2)
    {}

    inline Event(const IntEvent& coll, const unsigned long& nCC2) throw():
      dt(coll.getdt()),
      collCounter2(nCC2),
      type(INTERACTION),
      p2(coll.getParticle2ID())
    {
      if (coll.getType() == RECALCULATE) type = RECALCULATE;
    }

    inline Event(const GlobalEvent& coll) throw():
      dt(coll.getdt()),
      type(GLOBAL),
      p2(coll.getGlobalID())
    {
      if (coll.getType() == RECALCULATE) type = RECALCULATE;
    }

    inline Event(const LocalEvent& coll) throw():
      dt(coll.getdt()),
      type(LOCAL),
      p2(coll.getLocalID())
    {
      if (coll.getType() == RECALCULATE) type = RECALCULATE;
    }

    inline bool operator< (const Event& ip) const throw()
    { return dt < ip.dt; }

    inline bool operator> (const Event& ip) const throw()
    { return dt > ip.dt; }

    inline void stream(const double& ndt) throw() { dt -= ndt; }

    mutable double dt;
    unsigned long collCounter2;
    EEventType type;
    size_t p2;  
  };
}
