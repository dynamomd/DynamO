/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "neighbourlist.hpp"
#include "../extcode/threadpool.hpp"

class SThreadedNBList: public CSNeighbourList
{
public:
  SThreadedNBList(const XMLNode&, DYNAMO::SimData* const);

  SThreadedNBList(DYNAMO::SimData* const, CSSorter*, size_t threadCount);

  virtual void addEvents(const CParticle&);
  
  virtual void operator<<(const XMLNode&);

  virtual void fullUpdate(const CParticle& part);

  virtual void fullUpdate(const CParticle& p1, const CParticle& p2);

  void streamParticles(const CParticle& part, const size_t& id) const;

  void threadAddLocalEvent(const CParticle& part, 
			   const size_t& id);

  void threadAddIntEvent(const CParticle& part, 
			 const size_t& id);

  void spawnThreadAddLocalEvent(const CParticle& part, 
				const size_t& id);

  void spawnThreadAddIntEvents(const CParticle& part, const size_t& id);

  void addEvents2(const CParticle& part, const size_t& id) const;

protected:
  virtual void outputXML(xmlw::XmlStream&) const;

  void addEventsInit(const CParticle&);

  void addGlobal(const CParticle& p1, const smrtPlugPtr<CGlobal>& glob);

  ThreadPool _threadPool;

  boost::mutex _sorterLock;
  
};
