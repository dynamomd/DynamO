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
#include <dynamo/schedulers/sorters/FEL.hpp>
#include <magnet/exception.hpp>
#include <magnet/xmlwriter.hpp>
#include <vector>
#include <cmath>

namespace dynamo {
  template<class PEL>
  class CBTFEL: public FEL
  {
  public:
    virtual void init(const size_t N) 
    {
      clear();
      _streamFreq = _N = N;
      _CBT.resize(2 * N);
      _Leaf.resize(N + 1, std::numeric_limits<size_t>::max());
      _Min.resize(N + 1);
      _eventCount.resize(N, 0);
    }

    void clear()
    {
      _CBT.clear();
      _Leaf.clear();
      _Min.clear();
      _N = 0;
      _NP = 0;
      _pecTime = 0.0;
      _streamFreq = 0;
      _nUpdate = 0; 
      _activeID = std::numeric_limits<size_t>::max();
      _eventCount.clear();
    }

    inline void stream(const double dt)
    {    
      _pecTime += dt;
      ++_nUpdate;

      if (!(_nUpdate % _streamFreq))
	{
	  for (auto& pDat : _Min)
	    pDat.stream(_pecTime);
	  _pecTime = 0.0;
	}
    }

    inline void invalidate(const size_t ID) {
      flushChanges(ID);
      //Blank approximately half the events by clearing the PEL of the
      //particle.
      _Min[ID+1].clear();
      //Catch the others with lazy deletion.
      ++_eventCount[ID];
    }

    inline void pop() {
      flushChanges(_CBT[1]-1);
      _Min[_CBT[1]].pop();
    }

    inline bool empty() {
      flushChanges();
      return _CBT.empty() || _Min[_CBT[1]].empty();
    }


    virtual Event top() {
      //empty() causes a flush
      if (empty()) M_throw() << "Event queue is empty!";
      Event next_event = _Min[_CBT[1]].top();
      while ((next_event._source == INTERACTION) 
	     && (next_event._particle2eventcounter != _eventCount[next_event._particle2ID])) {
	pop();
	//empty() causes a flush
	if (empty()) M_throw() << "Event queue is empty!";
	next_event = _Min[_CBT[1]].top();
      }

      next_event._dt -= _pecTime;
      return next_event;
    }

    inline void push(Event event)
    {
      //Exit early
#ifdef DYNAMO_DEBUG
      if (std::isnan(event._dt))
	M_throw() << "NaN value pushed into the sorter.";
#endif
      //Only push events which will actually happen
      if (event._dt != HUGE_VAL) {
	flushChanges(event._particle1ID);
	event._dt += _pecTime;
	if (event._source == INTERACTION)
	  event._particle2eventcounter = _eventCount[event._particle2ID];
	_Min[event._particle1ID + 1].push(event);
      }
    }

    inline void rescaleTimes(const double factor)
    {
      for (auto& pDat : _Min)
	pDat.rescaleTimes(factor);
      _pecTime *= factor;
    }

    protected:
    size_t _activeID;

    virtual void flushChanges(const size_t ID = std::numeric_limits<size_t>::max()) {
      if ((_activeID != ID) && (_activeID !=std::numeric_limits<size_t>::max()))
	{
	  if (_Min[_activeID + 1].empty() || (_Min[_activeID + 1].top()._dt == HUGE_VAL)) {
	    if (_Leaf[_activeID + 1] != std::numeric_limits<size_t>::max()) {
	      Delete(_activeID + 1);
	    }
	  } else {
	    if (_Leaf[_activeID + 1] == std::numeric_limits<size_t>::max()) {
	      Insert(_activeID + 1);
	    }
	    else {
	      UpdateCBT(_activeID + 1);
	    }
	  }
	}
      _activeID = ID;
    }
  
    std::vector<size_t> _CBT;
    std::vector<size_t> _Leaf;
    std::vector<PEL> _Min;
    size_t _NP, _N, _streamFreq, _nUpdate;
  
    double _pecTime;
  
    std::vector<size_t> _eventCount;

    ///////////////////////////BINARY TREE IMPLEMENTATION
    inline void UpdateCBT(const size_t i)
    {
      size_t f = _Leaf[i] / 2;
    
      for(; (f > 0) && (_CBT[f] == i); f /= 2) 
	{
	  size_t l = _CBT[f*2],
	    r = _CBT[f*2+1];
	  _CBT[f] = (_Min[r] > _Min[l]) ? l : r;
	}

      //Walk up finding the winners till it doesn't change or you hit
      //the top of the tree
      for( ; f>0; f /= 2) 
	{
	  size_t w = _CBT[f], /* old winner */
	    l = _CBT[f*2],
	    r = _CBT[f*2+1];
	
	  _CBT[f] = (_Min[r] > _Min[l]) ? l : r;

	  if (_CBT[f] == w) return; /* end of the event time comparisons */
	}
    }


    inline void Insert(const size_t i)
    {
      if (_NP)
	{
	  size_t j = _CBT[_NP];
	  _CBT[_NP*2] = j;
	  _CBT[_NP*2+1] = i;
	  _Leaf[j] = _NP*2;
	  _Leaf[i]= _NP*2+1;
	  ++_NP;
	  UpdateCBT(j);
	}
      else
	{
	  _CBT[1]=i;
	  _Leaf[i]=1;
	  ++_NP; 
	}
    }
  
    inline void Delete(const size_t i)
    {
      if (_NP < 2) { _CBT[1]=0; _Leaf[0]=1; --_NP; return; }

      size_t l = _NP * 2 - 1;

      if (_CBT[l-1] == i)
	{
	  _Leaf[_CBT[l]] = l/2;
	  _CBT[l/2] = _CBT[l];
	  UpdateCBT(_CBT[l]);
	  --_NP;
	  _Leaf[i] = std::numeric_limits<size_t>::max();
	  return;
	}

      _Leaf[_CBT[l-1]] = l/2;
      _CBT[l/2] = _CBT[l-1];
      UpdateCBT(_CBT[l-1]);

      if (_CBT[l] != i) 
	{
	  _CBT[_Leaf[i]] = _CBT[l];
	  _Leaf[_CBT[l]] = _Leaf[i];
	  UpdateCBT(_CBT[l]);
	}
    
      --_NP;
      _Leaf[i] = std::numeric_limits<size_t>::max();
    }

    virtual void outputXML(magnet::xml::XmlStream& XML) const
    { XML << magnet::xml::attr("Type") << (std::string("CBT") + PEL::name()); }
    };
  }
