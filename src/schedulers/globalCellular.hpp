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

#ifndef CSGlobCellular_H
#define CSGlobCellular_H

#include "scheduler.hpp"
#include <algorithm>
#include <queue>
#include <boost/foreach.hpp>
#include "sorters/boundedPQ.hpp"
#include "sorters/cbt.hpp"

//Use to disable the bounded priority queue when debugging bad events
//#define CBT

class CSGlobCellular: public CScheduler
{
public:
  friend class COPBoundedQStats;

  CSGlobCellular(const XMLNode&, const DYNAMO::SimData*);

  CSGlobCellular(const DYNAMO::SimData*);

  virtual void rebuildList() { initialise(); }

  virtual void initialise();

  virtual void notifyVirtualCellsReinit();

  virtual void update(const CParticle&);

  virtual ENextEvent nextEventType() const;

  virtual const CLocalEvent earliestLocalEvent() const;

  virtual const CGlobEvent earliestGlobEvent() const;

  virtual const CIntEvent earliestIntEvent() const;

  virtual void stream(const Iflt);

  virtual void rescaleTimes(Iflt);

  virtual void popVirtualEvent();

  void virtualCellNewNeighbour(const CParticle&, const CParticle&); 

  virtual void pushAndUpdateVirtualEvent(const CParticle&, const intPart&);

  virtual void operator<<(const XMLNode&);

protected:
  virtual void outputXML(xmlw::XmlStream&) const;

  void addNewEvents(const CParticle&) const;

#ifndef CBT
  mutable CSSBoundedPQ eventHeap;
# else
  mutable CSSCBT  eventHeap;
#endif

  mutable std::vector<unsigned long long> eventCount;

  size_t GlobCellID;
};

#endif
