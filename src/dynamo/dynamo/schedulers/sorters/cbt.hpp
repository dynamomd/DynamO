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
#include <dynamo/schedulers/sorters/datastruct.hpp>
#include <dynamo/schedulers/sorters/sorter.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <magnet/exception.hpp>
#include <magnet/xmlwriter.hpp>
#include <vector>
#include <cmath>

namespace dynamo {
  class CSSCBT: public CSSorter
  {
  private:
    std::vector<unsigned long> CBT;
    std::vector<unsigned long> Leaf;
    std::vector<pList> Min;
    unsigned long NP, N, streamFreq, nUpdate;

    double pecTime;

  public:  
    CSSCBT(const dynamo::SimData* const& SD):
      CSSorter(SD, "CBT")
    {}

    typedef std::vector<pList>::iterator iterator;
    typedef std::vector<pList>::const_iterator const_iterator;

    inline iterator begin() { return Min.begin(); }
    inline const_iterator begin() const { return Min.begin(); }
    inline iterator end() { return Min.end(); }
    inline const_iterator end() const { return Min.end(); }
    inline size_t size() const { return Min.size(); }
    inline bool empty() const { return Min.empty(); }

    void resize(const size_t& a)
    {
      clear();
      streamFreq = N = a;
      CBT.resize(2 * N);
      Leaf.resize(N + 1);
      Min.resize(N + 1);
    }

    void clear() 
    {
      CBT.clear();
      Leaf.clear();
      Min.clear();
      N = 0;
      NP = 0;
      pecTime = 0.0;
      streamFreq = 0;
      nUpdate = 0;
    }

    void init()
    {
      for (unsigned long i = 1; i <= N; i++)
	Insert(i);  
    }

    void rebuild()
    {
      init();
    }

    inline void stream(const double& dt)
    {    
      pecTime += dt;
      ++nUpdate;

      if (!(nUpdate % streamFreq))
	{
#ifdef dynamo_UpdateCollDebug
	  std::cerr << "PecTime Stream occuring";
#endif
	  BOOST_FOREACH(pList& pDat, Min)
	    BOOST_FOREACH(intPart& event, pDat)
	    event.dt -= pecTime;
	
	  pecTime = 0.0;
	}
    }

    inline void clearPEL(const size_t& ID) { Min[ID+1].clear(); }
    inline void popNextPELEvent(const size_t& ID) { Min[ID+1].pop(); }
    inline void popNextEvent() { Min[CBT[1]].pop(); }
    inline bool nextPELEmpty() const { return Min[CBT[1]].empty(); }

    inline intPart copyNextEvent() const 
    { intPart retval(Min[CBT[1]].top());
      retval.dt -= pecTime;
      return retval; 
    }

    inline EEventType next_type() const { return Min[CBT[1]].top().type; }
    inline unsigned long next_collCounter2() const { return Min[CBT[1]].top().collCounter2; }
    inline size_t next_p2() const { return Min[CBT[1]].top().p2; }

    inline void push(const intPart& tmpVal, const size_t& pID)
    {
      //Exit early
#ifdef DYNAMO_DEBUG
      if (boost::math::isnan(tmpVal.dt))
	M_throw() << "NaN value pushed into the sorter! Should be Inf I guess?";
#endif

      if (tmpVal.type == NONE) return;
      tmpVal.dt += pecTime;
      Min[pID+1].push(tmpVal);
    }

    inline void update(const size_t& a) { UpdateCBT(a+1); }

    inline double next_dt() const { return Min[CBT[1]].getdt() - pecTime; }

    inline size_t next_ID() const { return CBT[1] - 1; }

    inline void rescaleTimes(const double& factor)
    {
      BOOST_FOREACH(pList& pDat, Min)
	BOOST_FOREACH(intPart& event, pDat)
        event.dt *= factor;

      pecTime *= factor;
    }

    inline void sort() {}

  private:
    inline void UpdateCBT(unsigned int i)
    {
      unsigned int f = Leaf[i]/2,l,r,w;

      //While i is at the top we must keep walking up, cause i could win
      //or could not
      for(; f > 0; f = f / 2) 
	{
	  if (CBT[f] != i) break; /* jumps to the next "for" */
	  l = CBT[f*2];
	  r = CBT[f*2+1];
	  if (Min[r] > Min[l] )
	    CBT[f] = l;
	  else
	    CBT[f] = r;
	}
      //Walk up finding the winners till it doesn't change or you hit
      //the top of the tree
      for( ; f>0; f=f/2) 
	{
	  w = CBT[f]; /* old winner */
	  l = CBT[f*2];
	  r = CBT[f*2+1];
	  if (Min[r] > Min[l])
	    CBT[f] = l;
	  else
	    CBT[f] = r;
	  if (CBT[f] == w ) return; /* end of the event time comparisons */
	}
    }

    inline void Insert(unsigned int i)
    {
      if (NP == 0) {CBT[1]=i; NP++; return;}
      int j = CBT[NP];
      CBT [NP*2] = j;
      CBT [NP*2+1] = i;
      Leaf[j] = NP*2;
      Leaf[i]= NP*2+1;
      ++NP;
      UpdateCBT (j);
    }
  
    inline void Delete(unsigned int i)
    {
      if (NP < 2) { CBT[1]=0; Leaf[0]=1; --NP; return; }
      int l = NP * 2 - 1;
      if (CBT[l-1] != i)
	{
	  Leaf[CBT[l-1]] = l/2;
	  CBT[l/2] = CBT[l-1];
	  UpdateCBT(CBT[l-1]);
	}
      else 
	{
	  Leaf[CBT[l]] = l/2;
	  CBT[l/2] =CBT[l];
	  UpdateCBT(CBT[l]);
	  --NP;
	  return;
	}

      if (CBT[l] != i) 
	{
	  CBT[Leaf[i]] = CBT[l];
	  Leaf[CBT[l]] = Leaf[i];
	  UpdateCBT(CBT[l]);
	}
    
      --NP;
    }

    virtual void outputXML(magnet::xml::XmlStream& XML) const
    { XML << magnet::xml::attr("Type") << "CBT"; }

  };
}
