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
#include <dynamo/schedulers/sorters/sorter.hpp>
#include <dynamo/schedulers/sorters/heapPEL.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/simulation.hpp>
#include <magnet/exception.hpp>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>

#ifdef DYNAMO_DEBUG 
#include <boost/math/special_functions/fpclassify.hpp>
#endif

namespace dynamo {
  template<size_t Size>
  class PELMinMax;

  class PELSingleEvent;

  template<class T> struct FELBoundedPQName;

  template<>
  struct FELBoundedPQName<PELHeap>
  {
    inline static std::string name() { return "BoundedPQ"; }
  };

  template<size_t I>
  struct FELBoundedPQName<PELMinMax<I> >
  {
    inline static std::string name() { return std::string("BoundedPQMinMax") + boost::lexical_cast<std::string>(I); }
  };

  template<>
  struct FELBoundedPQName<PELSingleEvent>
  {
    inline static std::string name() { return "BoundedPQSingleEvent"; }
  };

  template<typename T = PELHeap>
  class FELBoundedPQ: public FEL
  {
  private:
    //Bounded priority queue variables and types
    struct eventQEntry
    {
      T data;
      int next;
      int previous;
      int qIndex;
    };

    std::vector<int> linearLists;
    int currentIndex;

    double scale;
    double pecTime;
    double listWidth;
    int nlists;  

    //Binary tree variables
    std::vector<unsigned long> CBT;
    std::vector<unsigned long> Leaf;
    std::vector<eventQEntry> Min;
    size_t NP, N;
    size_t exceptionCount;

  public:  
    FELBoundedPQ():exceptionCount(0) {}

    ~FELBoundedPQ()
    { 
      std::cout << "Exception Events = " << exceptionCount << std::endl;
    }

    void resize(const size_t& a)
    {
      clear();
      N = a;
      CBT.resize(2 * N);
      Leaf.resize(N + 1);
      Min.resize(N + 1);
      //Min.front().data.push(Event(HUGE_VAL, NONE));
    }

    void clear() 
    {
      Min.clear();
      CBT.clear();
      Leaf.clear();
      linearLists.clear();
      N = 0;
      NP = 0;
      currentIndex = 0;
      pecTime = 0.0;
    }

    inline void stream(const double& ndt) { pecTime += ndt; }

    void init()
    {
      init(false);
    }

    void rebuild()
    {
      init(true);
    }

    void init(bool quiet)
    {
      NP = 0;
      
      {
	double minVal(0), maxVal(-HUGE_VAL);
	size_t counter(0);
	
	try {
	  //Determine nlists and scale by instrumenting the queue
	  //We make a bounded queue of length N
	  //With a date length of the mean time between events
	  
	  for (const eventQEntry& dat : Min)
	    if (!std::isinf(dat.data.getdt()))
	      {
		minVal = std::min(minVal, dat.data.getdt());
		maxVal = std::max(maxVal, dat.data.getdt());
		++counter;
	      }
	  
	  if (counter < 10)
	    {
	      //Something is peculiar about the system
	      std::cerr <<
		"The event queue doesn't have more than 10 VALID events in it"
		"\nThis means the queue cannot be instrumented properly to"
		"\ndetermine the optimal settings for the bounded queue, now"
		"\nusing some (probably inefficient) defaults."
		"\nIf this is a proper simulation, consider using a different Sorter (e.g., CBT)." 
		   << std::endl;
	      scale = 10;
	      nlists = 1000;
	    }
	  else
	    {
	      if (maxVal < 0 ) 
		std::cerr << "WARNING! The event queue is filled with negative events!"  
		     << std::endl;;
	      
	      scale = counter / (maxVal - minVal);
	      nlists = Min.size();
	    }
	}
	catch (std::exception& excep)
	  {
	    M_throw() << "Failure in boundedPQ init()"
		      << "\nminVal = " << minVal
		      << "\nmaxVal = " << maxVal
		      << "\ncounter = " << counter
		      << "\n" << excep.what();
	  }
      }

      listWidth = nlists / scale;
      if (scale == HUGE_VAL)
	M_throw() << "The scale factor for the bounded priority queue is infinite. Cannot resolve this. May be caused by only having zero time collisions.";

      if (scale <= 0.0)
	M_throw() << "The scale factor for the bounded priority queue is zero. Cannot resolve this. May be caused by a large number of negative time events.";

      if (nlists == 0)
	{
	  std::cerr << "nlists = 0!\n"
		    << "This is a BAD thing, unless NCells = NParticles and "
	    "they're in a perfect crystal, if it happens again after the "
	    "preliminary run its certainly a bug" << std::endl;
	  nlists = 1000;
	}

      if (!quiet)
	std::cout << "Length of linear list = " << nlists
		  << "Scale factor = " << scale
		  << std::endl;

      linearLists.clear();
      linearLists.resize(nlists+1, -1); /*+1 for overflow, -1 for
					  marking empty*/ 

      if (!quiet)
	std::cout << "Sorting all events, please wait..." << std::endl;

      //Now insert all of the events!
      for (unsigned long i = 1; i <= N; i++)
	insertInEventQ(i);

  
      if (!quiet)
	std::cout << "Finding first event..." << std::endl;
    
      //Find the next event and place it first so nextEventID() works
      orderNextEvent();
      if (!quiet)
	std::cout << "Ready for simulation." << std::endl;
    }

    inline void push(const Event& tmpVal, const size_t& pID)
    {
#ifdef DYNAMO_DEBUG
      if (boost::math::isnan(tmpVal.dt))
	M_throw() << "NaN value pushed into the sorter! Should be Inf I guess?";
#endif 

      tmpVal.dt += pecTime;
      Min[pID + 1].data.push(tmpVal);
    }

    inline void update(const size_t& pID)
    {
      deleteFromEventQ(pID + 1);
      insertInEventQ(pID + 1);
    }

    inline void clearPEL(const size_t& ID) { Min[ID+1].data.clear(); }
    inline void popNextPELEvent(const size_t& ID) { Min[ID+1].data.pop(); }
    inline void popNextEvent() { Min[CBT[1]].data.pop(); }
    inline bool nextPELEmpty() const { return Min[CBT[1]].data.empty(); }

    inline size_t next_ID() const { return CBT[1] - 1; }
    inline EEventType next_type() const { return Min[CBT[1]].data.top().type; }
    inline unsigned long next_collCounter2() const { return Min[CBT[1]].data.top().collCounter2; }
    inline size_t next_p2() const { return Min[CBT[1]].data.top().p2; }

    inline double next_dt() const { return Min[CBT[1]].data.getdt() - pecTime; }

    inline void sort() { orderNextEvent(); }

    inline void rescaleTimes(const double& factor)
    {
      for (eventQEntry& dat : Min)
	dat.data.rescaleTimes(factor);

      pecTime *= factor;
    
      scale /= factor;
      listWidth = nlists/scale;

    }

  private:
    ///////////////////////////BOUNDED QUEUE IMPLEMENTATION
    inline void insertInEventQ(int p)
    {
      double box = scale * Min[p].data.getdt();

      int i = (box > std::numeric_limits<int>::max())
	? (nlists + nlists) //Put this in the overflow list
	: static_cast<int>(box); //You can use this as usual
    
      //This line makes negative time events possible without a segfault
      i = std::max(i, currentIndex);

      if (i > (nlists-1)) /* account for wrap */
	{
	  i -= nlists;
	  if(i>=currentIndex-1)
	    //Its overflowed!
	    i=nlists; /* store in overflow list */
	}

      Min[p].qIndex=i;

      if(i == currentIndex)
	Insert(p); /* insert in PQ */
      else
	{
	  /* insert in linked list */
	  int oldFirst = linearLists[i];
	  Min[p].previous = -1;
	  Min[p].next = oldFirst;
	  linearLists[i]= p;
	  if(oldFirst != -1)
	    Min[oldFirst].previous = p;
	}
    }

    inline void processOverflowList()
    {
      int e = linearLists[nlists];
      linearLists[nlists] = -1; /* mark empty; we will treat all entries and may re-add some */

      size_t overflowEvents = 0;
      while(e!=-1)
	{
	  ++overflowEvents;
	  int eNext = Min[e].next; /* save next */
	  insertInEventQ(e); /* try add to regular list now */
	  e=eNext;
	}

      exceptionCount += overflowEvents;
      
      //Check if the overflow contained more than half the total
      //events. If so, force a complete rebuild of the scheduler
      if (overflowEvents > N/2)
	{
	  init();
	}
    }

    inline void deleteFromEventQ(const int& e)
    {
      if(Min[e].qIndex == currentIndex)
	Delete(e); /* delete from pq */
      else
	{
	  /* remove from linked list */
	  int prev = Min[e].previous,
	    next = Min[e].next;
	  if(prev == -1)
	    linearLists[Min[e].qIndex] = Min[e].next;
	  else
	    Min[prev].next = next;

	  if(next != -1)
	    Min[next].previous = prev;
	}
    }

    inline void orderNextEvent()
    {
      while(NP==0)
	{
	  /*The current priority queue is exhausted, move on to the
	    next one*/

	  /* change current calendar "date" */
	  if(++currentIndex==nlists)
	    {
	      /* We've reached the last "date" in the calendar.
	       Reset the index (wrap the date).*/
	      currentIndex=0;

	      //Stream every event by the list width!
	      for (eventQEntry& dat : Min)
		dat.data.stream(listWidth);
	      //update the peculiar time
	      pecTime -= listWidth;

	      //Need to process this once per wrap so do it now 
	      //All events that had dt > listWidth are now processed
	      processOverflowList();
	    }

	  /* populate pq */
	  for (int e = linearLists[currentIndex]; e!=-1; e=Min[e].next)
	    Insert(e);

	  linearLists[currentIndex] = -1;
	}
    }


    ///////////////////////////BINARY TREE IMPLEMENTATION
    inline void UpdateCBT(const unsigned int& i)
    {
      unsigned int f = Leaf[i] / 2;
    
      for(; (f > 0) && (CBT[f] == i); f /= 2) 
	{
	  unsigned int l = CBT[f*2],
	    r = CBT[f*2+1];
	  CBT[f] = (Min[r].data > Min[l].data) ? l : r;
	}

      //Walk up finding the winners till it doesn't change or you hit
      //the top of the tree
      for( ; f>0; f /= 2) 
	{
	  unsigned int w = CBT[f], /* old winner */
	    l = CBT[f*2],
	    r = CBT[f*2+1];
	
	  CBT[f] = (Min[r].data > Min[l].data) ? l : r;

	  if (CBT[f] == w) return; /* end of the event time comparisons */
	}
    }

    inline void Insert(const unsigned int& i)
    {
      if (NP)
	{
	  int j = CBT[NP];
	  CBT [NP*2] = j;
	  CBT [NP*2+1] = i;
	  Leaf[j] = NP*2;
	  Leaf[i]= NP*2+1;
	  ++NP;
	  UpdateCBT(j);
	}
      else
	{
	  CBT[1]=i;
	  ++NP; 
	}
    }
  
    inline void Delete(const unsigned int& i)
    {
      if (NP < 2) { CBT[1]=0; Leaf[0]=1; --NP; return; }

      int l = NP * 2 - 1;

      if (CBT[l-1] == i)
	{
	  Leaf[CBT[l]] = l/2;
	  CBT[l/2] =CBT[l];
	  UpdateCBT(CBT[l]);
	  --NP;
	  return;
	}

      Leaf[CBT[l-1]] = l/2;
      CBT[l/2] = CBT[l-1];
      UpdateCBT(CBT[l-1]);

      if (CBT[l] != i) 
	{
	  CBT[Leaf[i]] = CBT[l];
	  Leaf[CBT[l]] = Leaf[i];
	  UpdateCBT(CBT[l]);
	}
    
      --NP;
    }
  
    virtual void outputXML(magnet::xml::XmlStream& XML) const
    { XML << magnet::xml::attr("Type") << FELBoundedPQName<T>::name(); }

  };
}
