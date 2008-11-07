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

#ifndef CSNeighbourList_H
#define CSNeighbourList_H

#include "scheduler.hpp"
#include <algorithm>
#include <queue>
#include <boost/foreach.hpp>
#include <boost/signals.hpp>

class CSNeighbourList: public CScheduler
{
public:
  friend class COPBoundedQStats;

  CSNeighbourList(const XMLNode&, const DYNAMO::SimData*);

  CSNeighbourList(const DYNAMO::SimData*);

  virtual void rebuildList() { initialise(); }

  virtual void initialise();

  virtual void update(const CParticle&);

  virtual ENextEvent nextEventType() const;

  virtual const CLocalEvent earliestLocalEvent() const;

  virtual const CGlobEvent earliestGlobEvent() const;

  virtual const CIntEvent earliestIntEvent() const;

  virtual void stream(const Iflt&);

  virtual void rescaleTimes(const Iflt&);

  virtual void popVirtualEvent();

  void virtualCellNewNeighbour(const CParticle&, const CParticle&); 

  virtual void pushAndUpdateVirtualEvent(const CParticle&, const intPart&);

  virtual void operator<<(const XMLNode&);

protected:
  virtual void outputXML(xmlw::XmlStream&) const;

  void addNewEvents(const CParticle&) const;
  
  mutable std::vector<unsigned long long> eventCount;

  void addInteractionEvent(const CParticle&, const size_t&) const;
  void addLocalEvent(const CParticle&, const size_t&) const;
  
  size_t NBListID;

  boost::signals::scoped_connection cellChange;
  boost::signals::scoped_connection reinit;  
};

#endif
