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
#include <dynamo/schedulers/sorters/CBTFEL.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/simulation.hpp>
#include <magnet/exception.hpp>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>

static const size_t NO_LINK = std::numeric_limits<size_t>::max();

namespace dynamo {
  namespace detail {
    template<class PEL>
    struct BPQEntry : public PEL {
      BPQEntry(): next(NO_LINK), previous(NO_LINK), qIndex(NO_LINK) {}
      size_t next, previous, qIndex;
    };
  }

  template<typename PEL>
  class BoundedPQFEL: public CBTFEL<detail::BPQEntry<PEL> >
  {
    typedef CBTFEL<detail::BPQEntry<PEL> > Base;
  private:
    //Bounded priority queue variables and types

    std::vector<size_t> linearLists;
    size_t currentIndex;

    double scale;
    size_t nlists;
    size_t exceptionCount;
    size_t _optimizeCounter;
    
  public:  
    BoundedPQFEL():exceptionCount(0) {}
    
    ~BoundedPQFEL() { 
      std::cout << "Exception Events = " << exceptionCount << std::endl;
    }

    void init(const size_t N)
    {
      clear();
      Base::init(N);
      Base::_streamFreq = 100;

      //Start with the FEL in CBT mode
      scale=0;
      nlists = 1;
      linearLists.resize(nlists+1, NO_LINK); /*+1 for overflow, NO_LINK for marking empty*/ 
    }

    void clear()
    {
      Base::clear();
      linearLists.clear();
      currentIndex = 0;
      _optimizeCounter = 1;
    }

    inline void stream(const double ndt) {
      Base::_pecTime += ndt; 
    }

    inline void rescaleTimes(const double factor)
    {
      for (auto& dat : Base::_Min)
	dat.rescaleTimes(factor);

      Base::_pecTime *= factor;
      scale /= factor;
    }

  private: 
    virtual void flushChanges(const size_t ID = std::numeric_limits<size_t>::max()) {
      if ((Base::_activeID != ID) && (Base::_activeID !=std::numeric_limits<size_t>::max()))
	{
	  ////Optimise the queue settings every 10^6 events or so
	  //if (!(++_optimizeCounter % 2^20)) {
	  //  optimiseSettings();
	  //  return;
	  //}

	  insertInEventQ(Base::_activeID + 1);
	  orderNextEvent();
	}
      Base::_activeID = ID;
    }

    void optimiseSettings() {
      //Collect statistics on the event list.
      double minVal(HUGE_VAL), maxVal(-HUGE_VAL);
      size_t counter(0);
      
      for (const auto& dat : Base::_Min)
	if (!std::isinf(dat.top()._dt))
	  {
	    minVal = std::min(minVal, dat.top()._dt);
	    maxVal = std::max(maxVal, dat.top()._dt);
	    ++counter;
	  }
      
      if ((counter < 10) || (maxVal < 0) || ((maxVal - minVal) == 0) || (maxVal < minVal))
	{ //In unusual systems, drop down to a CBT queue
	  scale = 0;
	  nlists = 1;
	}
      else
	{
	  scale = counter / (maxVal - minVal);
	  nlists = Base::_Min.size();
	}

      //Mark all PELs as uninserted
      Base::_NP = 0;
      linearLists.clear();
      linearLists.resize(nlists+1, NO_LINK); /*+1 for overflow, NO_LINK for marking empty*/ 

      //Now insert all PELs
      for (unsigned long i = 1; i <= Base::_N; i++) {
	Base::_Min[i].qIndex = NO_LINK;
	Base::_Min[i].next = NO_LINK;
	Base::_Min[i].previous = NO_LINK;
	insertInEventQ(i);
      }

      orderNextEvent();
    }


    ///////////////////////////BOUNDED QUEUE IMPLEMENTATION
    inline void insertInEventQ(const size_t p)
    {
      //If its already inserted, then delete it first
      if (Base::_Min[p].qIndex != NO_LINK)
	deleteFromEventQ(p);

      //Check that the Q is not empty or filled with events which will never happen
      if (Base::_Min[p].empty() || (Base::_Min[p].top()._dt == HUGE_VAL))
	//Don't bother adding it to the queue.
	return;

      const double box = scale * Base::_Min[p].top()._dt;

      size_t i = (box > std::numeric_limits<size_t>::max())
	? (nlists + nlists) //Put this in the overflow list
	: static_cast<size_t>(box); //You can use this as usual
    
      //This line makes negative time events possible without a segfault
      i = std::max(i, currentIndex);

      if (i > (nlists-1)) /* account for wrap */
	{
	  i -= nlists;
	  if(i>=currentIndex-1)
	    //Its overflowed!
	    i=nlists; /* store in overflow list */
	}

      Base::_Min[p].qIndex=i;

      if(i == currentIndex)
	Base::Insert(p); /* insert in PQ */
      else
	{
	  /* insert in linked list */
	  size_t oldFirst = linearLists[i];
	  Base::_Min[p].previous = NO_LINK;
	  Base::_Min[p].next = oldFirst;
	  linearLists[i]= p;
	  if(oldFirst != NO_LINK)
	    Base::_Min[oldFirst].previous = p;
	}
    }

    inline void processOverflowList()
    {
      size_t e = linearLists[nlists];
      linearLists[nlists] = NO_LINK; /* mark empty; we will treat all entries and may re-add some */

      size_t overflowEvents = 0;
      while(e!=NO_LINK)
	{
	  ++overflowEvents;
	  size_t eNext = Base::_Min[e].next; /* save next */
	  insertInEventQ(e); /* try add to regular list now */
	  e = eNext;
	}
      exceptionCount += overflowEvents;
    }

    inline void deleteFromEventQ(const size_t e)
    {
      if(Base::_Min[e].qIndex == currentIndex)
	Base::Delete(e); /* delete from pq */
      else if (Base::_Min[e].qIndex != NO_LINK) {
	/* remove from linked list */
	const size_t prev = Base::_Min[e].previous,
	  next = Base::_Min[e].next;
	if(prev == NO_LINK)
	  linearLists[Base::_Min[e].qIndex] = Base::_Min[e].next;
	else
	  Base::_Min[prev].next = next;
	
	if(next != NO_LINK)
	  Base::_Min[next].previous = prev;
      }
      
      Base::_Min[e].qIndex = NO_LINK;
    }

    inline void orderNextEvent()
    {
      while(Base::_NP==0)
	{
	  /*The current priority queue is exhausted, move on to the
	    next one*/

	  /* change current calendar "date" */
	  if(++currentIndex == nlists)
	    {
	      /* We've reached the last "date" in the calendar.
	       Reset the index (wrap the date).*/
	      currentIndex = 0;

	      //Stream every event by the list width, and check if the
	      //queue is actually empty.
	      bool no_events = true;
	      const double listWidth = nlists / scale;
	      for (auto& dat : Base::_Min) {
		no_events = no_events && ((dat.empty()) || (dat.top()._dt == HUGE_VAL));
		dat.stream(listWidth);
	      }
	      //update the peculiar time
	      Base::_pecTime -= listWidth;
		
	      //Check if there are no events to schedule!
	      if (no_events && (linearLists[nlists] == NO_LINK))
		return;

	      //Need to process this once per wrap so do it now 
	      //All events that had dt > listWidth are now processed
	      processOverflowList();
	    }

	  /* populate pq */
	  for (size_t e = linearLists[currentIndex]; e != NO_LINK; e = Base::_Min[e].next)
	    Base::Insert(e);
	  linearLists[currentIndex] = NO_LINK;
	}
    }

  
    virtual void outputXML(magnet::xml::XmlStream& XML) const { 
      XML << magnet::xml::attr("Type") << (std::string("BoundedPQ") + PEL::name()); 
    }

  };
}
