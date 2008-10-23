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

#ifndef SCHEDULER_H
#define SCHEDULER_H

#include "../base/is_base.hpp"
#include "sorters/cbt.hpp"
#include "../dynamics/interactions/intEvent.hpp"
#include "../dynamics/globals/globEvent.hpp"
#include <vector>

typedef enum {
  Interaction , /*!< No collision occurs*/
  Global      , /*!< Hard core collision*/
  System        /*!< Well Event, could be Dissociation/bounce etc*/
} ENextEvent;

class CParticle;

class CScheduler: public DYNAMO::SimBase_const
{
public:
  CScheduler(const DYNAMO::SimData* const tmp, const char * aName):
    SimBase_const(tmp, aName, IC_purple)
  {}
  
  virtual ~CScheduler() = 0;

  virtual void initialise() = 0;

  virtual void update(const CParticle&) = 0;

  inline void update(const CParticle& part1, const CParticle& part2)
  {
    update(part1);
    update(part2);
  }

  virtual void popVirtualEvent() = 0;

  virtual void virtualCellNewNeighbour(const CParticle&, const CParticle&) =0; 

  virtual void pushAndUpdateVirtualEvent(const CParticle&, const intPart&) = 0;

  virtual void stream(const Iflt) = 0;
  
  virtual const CIntEvent earliestIntEvent() const = 0;

  virtual const CGlobEvent earliestGlobEvent() const = 0;

  virtual ENextEvent nextEventType() const = 0;

  virtual void rebuildList() = 0;

  friend xmlw::XmlStream& operator<<(xmlw::XmlStream&, const CScheduler&);

  static CScheduler* getClass(const XMLNode&, const DYNAMO::SimData*);

  virtual void operator<<(const XMLNode&) = 0;
  
  virtual void rescaleTimes(Iflt) { D_throw() << "Not implemented yet"; }

protected:

  virtual void outputXML(xmlw::XmlStream&) const = 0;

  CGlobEvent getGlobEvent(const CParticle&) const;
};

#endif
