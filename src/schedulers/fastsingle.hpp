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

#ifndef CSFastSingle_H
#define CSFastSingle_H

#include "scheduler.hpp"
#include <vector>
#include "sorters/generalcbt.hpp"

class CSFastSingle : public CScheduler
{
 public:
  CSFastSingle(const XMLNode&, const DYNAMO::SimData*);
  CSFastSingle(const DYNAMO::SimData*);

  virtual void initialise();

  virtual void update(const CParticle&);

  virtual ENextEvent nextEventType() const;

  virtual void rebuildList();

  virtual void operator<<(const XMLNode&);

  virtual void rescaleTimes(Iflt);

  const CIntEvent earliestIntEvent() const;
  const CGlobEvent earliestGlobEvent() const;
  const CLocalEvent earliestLocalEvent() const
  { D_throw() << "Not implemented"; }

  void stream(const Iflt);

  virtual void popVirtualEvent();

  virtual void pushAndUpdateVirtualEvent(const CParticle&, const intPart&);

 protected:
  virtual void outputXML(xmlw::XmlStream&) const;
  void initGlobalQueue();

  mutable std::vector<CIntEvent> intEventQueue;
  mutable std::vector<CGlobEvent> globEventQueue;
  mutable std::vector<CIntEvent>::const_iterator nextIntEvent;
  mutable std::vector<CGlobEvent>::const_iterator nextGlobEvent; 

 private:
  void updateCollision(CIntEvent &);
  void rebuildCollision(CIntEvent &);
};

#endif
