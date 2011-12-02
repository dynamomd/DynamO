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
#include <dynamo/schedulers/sorters/sorter.hpp>
#include <dynamo/schedulers/sorters/heapPEL.hpp>
#include <dynamo/dynamics/units/units.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <boost/static_assert.hpp>
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

  template<class T>
  struct FELBoundedPQName
  {
    BOOST_STATIC_ASSERT(sizeof(T) == 0);
  };

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
    FELBoundedPQ(const dynamo::SimData* const& SD):
      FEL(SD, "BoundedPQ"),
      exceptionCount(0) 
    {}

    ~FELBoundedPQ() 
    { 
      dout << "Exception Events = " << exceptionCount << std::endl;
    }
  
    inline size_t size() const { return Min.size() - 1; }
    inline bool empty() const { return Min.empty(); }

    inline const int& NLists() const { return nlists; }
    inline const double& scaleFactor() const { return scale; }
    inline const size_t& exceptionEvents() const { return exceptionCount; }
    inline const size_t& treeSize() const { return NP; }

    inline std::vector<size_t> getEventCounts() const
    {
      std::vector<size_t> tmpVec;
      tmpVec.resize(nlists - 1,0);

      //Miss the binary tree
      for (int i = 1; i < nlists - 1; ++i)
	{
	  int index = i + currentIndex;
	  if (index > nlists -1)
	    index -= nlists;
	  size_t counter = 0;
	  //Scroll through the list counting
	  int nextID = linearLists[index];	
	  while (nextID != -1)
	    {
	      ++counter;
	      nextID = Min[nextID].next;
	    }
	  tmpVec[i] = counter;
	}
      return tmpVec;
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
      CBT.clear();
      Leaf.clear();
      Min.clear();
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
      {
	double minVal(0), maxVal(-HUGE_VAL);
	size_t counter(0);
	
	try {
	  //Determine nlists and scale by instrumenting the queue
	  //We make a bounded queue of length N
	  //With a date length of the mean time between events
	  
	  BOOST_FOREACH(const eventQEntry& dat, Min)
	    if (!std::isinf(dat.data.getdt()))
	      {
		if (dat.data.getdt() < minVal) minVal = dat.data.getdt();
		if (dat.data.getdt() > maxVal) maxVal = dat.data.getdt();
		++counter;
	      }
	  
	  if (counter < 10)
	    {
	      //Something is peculiar about the system
	      derr <<
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
		derr << "WARNING! The event queue is filled with negative events!"  
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
	M_throw() << "Scale factor is infinite (only zero time collisions or no collisions?)";

      if (scale <= 0.0)
	M_throw() << "Scale factor is zero or negative (negative collisions?)";

      if (nlists == 0)
	{
	  derr << "nlists = 0!\n"
	       << "This is a BAD thing, unless NCells = NParticles and "
	    "they're in a perfect crystal, if it happens again after the "
	    "preliminary run its certainly a bug" << std::endl;
	  nlists = 1000;
	}

      if (!quiet)
	dout << "Length of linear list = " << nlists
	     << "Scale factor = " 
	     << scale * Sim->dynamics.units().unitTime() 
	     << std::endl;

      linearLists.resize(nlists+1, -1); /*+1 for overflow, -1 for
					  marking empty*/ 

      if (!quiet)
	{
	  dout << "Sorting all events, please wait..." << std::endl;
	}

      //Now insert all of the events!
      for (unsigned long i = 1; i <= N; i++)
	insertInEventQ(i);

  
      if (!quiet)
	{
	  dout << "Finding first event..." << std::endl;
	}
    
      //Find the next event and place it first so nextEventID() works
      orderNextEvent();
      if (!quiet)
	dout << "Ready for simulation." << std::endl;
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

    //  inline const T& operator[](const size_t& a) const 
    //  {
    //#ifdef DYNAMO_DEBUG 
    //    if (Min.empty())
    //      M_throw() << "Heap not yet sized";
    //#endif
    //    
    //    return Min[a+1].data; 
    //  }
  
    //  inline T& operator[](const size_t& a) 
    //  {
    //#ifdef DYNAMO_DEBUG 
    //    if (Min.empty())
    //      M_throw() << "Heap not yet sized";
    //#endif
    //
    //    return Min[a+1].data; 
    //  }

    inline void clearPEL(const size_t& ID) { Min[ID+1].data.clear(); }
    inline void popNextPELEvent(const size_t& ID) { Min[ID+1].data.pop(); }
    inline void popNextEvent() { Min[CBT[1]].data.pop(); }
    inline bool nextPELEmpty() const { return Min[CBT[1]].data.empty(); }

    inline Event copyNextEvent() const 
    { Event retval(Min[CBT[1]].data.top());
      retval.dt -= pecTime;
      return retval; 
    }

    inline size_t next_ID() const { return CBT[1] - 1; }
    inline EEventType next_type() const { return Min[CBT[1]].data.top().type; }
    inline unsigned long next_collCounter2() const { return Min[CBT[1]].data.top().collCounter2; }
    inline size_t next_p2() const { return Min[CBT[1]].data.top().p2; }

    //inline T& next_Data() { return Min[CBT[1]].data; }
    //inline const T& next_Data() const { return Min[CBT[1]].data; }
    inline double next_dt() const { return Min[CBT[1]].data.getdt() - pecTime; }

    inline void sort() { orderNextEvent(); }

    inline void rescaleTimes(const double& factor)
    {
      BOOST_FOREACH(eventQEntry& dat, Min)
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
      if (i < currentIndex) i = currentIndex;

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

      while(e!=-1)
	{
	  ++exceptionCount;
	  int eNext = Min[e].next; /* save next */
	  insertInEventQ(e); /* try add to regular list now */
	  e=eNext;
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
      while(NP==0)/*if priority queue exhausted*/
	{
#ifdef dynamo_UpdateCollDebug
	  std::cerr << "\nQueue exhausted";
#endif
	  /* change current index */
	  if(++currentIndex==nlists)
	    {
	      //This is where we've wrapped all the way around
	      //Need to do lots here to maintain the list
	      currentIndex=0;

	      //Stream every event by the list width!
	      BOOST_FOREACH(eventQEntry& dat, Min)
		dat.data.stream(listWidth);
#ifdef dynamo_UpdateCollDebug
	      std::cerr << "\nPecTime Stream occuring";
#endif
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
