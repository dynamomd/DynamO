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

#ifndef CSSBoundedPQ_H
#define CSSBoundedPQ_H

#include <vector>
#include <cmath>
#include <iostream>
#include "../../base/is_exception.hpp"
#include "datastruct.hpp"
#include "sorter.hpp"

class CSSBoundedPQ: public CSSorter
{
private:
  //Bounded priority queue variables and types
  struct eventQEntry
  {
    int next;
    int previous;
    pList data;
    int qIndex;
  };

  std::vector<int> linearLists;
  int currentIndex;

  Iflt scale;
  Iflt pecTime;
  Iflt listWidth;
  int nlists;  

  //Binary tree variables
  std::vector<unsigned long> CBT;
  std::vector<unsigned long> Leaf;
  std::vector<eventQEntry> Min;
  size_t NP, N;
  size_t exceptionCount;

public:  
  CSSBoundedPQ():exceptionCount(0) {}
  ~CSSBoundedPQ() 
  { std::cout << "\nBPQ: Exception Events = " << exceptionCount << "\n"; }
  
  inline size_t size() const { return Min.size() - 1; }
  inline bool empty() const { return Min.empty(); }

  inline size_t NLists() { return nlists; }
  inline Iflt scaleFactor() { return scale; }
  inline Iflt exceptionEvents() { return exceptionCount; }
  inline size_t treeSize() { return NP; }

  inline std::vector<size_t> getEventCounts()
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
    Min.front().data.push(intPart(HUGE_VAL, NONE));
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

  inline void stream(const Iflt& ndt) { pecTime += ndt; }

  void init()
  {
    //Determine nlists and scale by instrumenting the queue
    std::vector<pList> tmpList;
    tmpList.reserve(Min.size());
    
    BOOST_FOREACH(eventQEntry& dat, Min)
      tmpList.push_back(dat.data);

    std::sort(tmpList.begin(), tmpList.end());

    //Find the mean dt value for events that are not infinite
    Iflt acc = 0.0;
    long counter = 0;

    for (std::vector<pList>::iterator iPtr = tmpList.begin() +1; iPtr != tmpList.end(); ++iPtr)
      if (iPtr->getdt() != HUGE_VAL)
	{ acc += iPtr->getdt() - (iPtr -1)->getdt(); counter++; }
      else
	break; //Jump out, nothing useful left

    if (counter < 2)
      {
	//Something is peculiar about the system
	std::cerr << IC_red << 
	  "BOUNDEDPQ: The event queue doesn't have more than 2 events in it"
	  "\nBOUNDEDPQ: This means the queue cannot be instrumented to"
	  "\nBOUNDEDPQ: determine the settings for the bounded queue, just"
	  "\nBOUNDEDPQ: using something that hopes the events in sim time"
	  "\nBOUNDEDPQ: arent longer than t=10000\n"
		  << IC_reset;

	init(10, 1000);
      }
    else
      {
	Iflt nscale = counter / acc;
	
	//Determine where the queue ends
	std::vector<pList>::reverse_iterator rIt = tmpList.rbegin();
	while (rIt->getdt() == HUGE_VAL)
	  ++rIt;
	
	//Determine the number of lists as the number required to cover
	//the current list plus some factor to reduce exceptions and
	//stream events
	int newnlists = static_cast<int>(2.0 * rIt->getdt() * nscale);
	
	init(nscale, newnlists);
      }
  }

  void init(const Iflt tmpScale, const int tmpnlists)
  {
    scale = tmpScale; 
    nlists = tmpnlists;
    listWidth = nlists / scale;
    if (scale == HUGE_VAL)
      D_throw() << "Scale factor is infinite (only zero time collisions or no collisions?)";

    if (scale <= 0.0)
      D_throw() << "Scale factor is zero or negative (negative collisions?)";

    if (nlists == 0)
      {
	std::cout << "\nBOUNDEDPQ: nlists = 0!";
	std::cout << "\nBOUNDEDPQ: This is a BAD thing, unless NCells = NParticles and they're in a perfect crystal, if it happens again after the preliminary run its a bug";
	nlists = 1000;
      }

    linearLists.resize(nlists+1, -1); /*+1 for overflow, -1 for
					marking empty*/ 

    //Now insert all of the events!
    for (unsigned long i = 1; i <= N; i++)
      insertInEventQ(i);

    //Find the next event and place it first so nextEventID() works
    orderNextEvent();
  }

  inline void push(const intPart& tmpVal, const size_t& pID)
  {
    tmpVal.dt += pecTime;
    Min[pID + 1].data.push(tmpVal);
  }

  inline void update(const int& pID)
  {
    deleteFromEventQ(pID + 1);
    insertInEventQ(pID + 1);
  }

  inline const pList& operator[](const size_t& a) const 
  {
#ifdef DYNAMO_DEBUG 
    if (Min.empty())
      D_throw() << "Heap not yet sized";
#endif
    
    return Min[a+1].data; 
  }
  
  inline pList& operator[](const size_t& a) 
  {
#ifdef DYNAMO_DEBUG 
    if (Min.empty())
      D_throw() << "Heap not yet sized";
#endif

    return Min[a+1].data; 
  }

  inline size_t next_ID() const { return CBT[1] - 1; }
  inline pList& next_Data() { return Min[CBT[1]].data; }
  inline const pList& next_Data() const { return Min[CBT[1]].data; }
  inline Iflt next_dt() const { return Min[CBT[1]].data.getdt() - pecTime; }

  inline void sort() { orderNextEvent(); }

  inline void rescaleTimes(const Iflt& factor)
  {
    BOOST_FOREACH(eventQEntry& dat, Min)
      dat.data.rescaleTimes(factor);

    pecTime *= factor;
    
    scale /= factor;
    listWidth = nlists/scale;

  }

private:
  virtual CSSorter* Clone() const { return new CSSBoundedPQ(*this); };
  ///////////////////////////BOUNDED QUEUE IMPLEMENTATION
  inline void insertInEventQ(int p)
  {
    Iflt box = scale * Min[p].data.getdt();

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
    int e = linearLists[nlists],
      eNext;
    linearLists[nlists] = -1; /* mark empty; we will treat all entries and may re-add some */

    while(e!=-1)
      {
	++exceptionCount;
	eNext=Min[e].next; /* save next */
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
#ifdef DYNAMO_UpdateCollDebug
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
#ifdef DYNAMO_UpdateCollDebug
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
};
#endif
