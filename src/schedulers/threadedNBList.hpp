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
#include <magnet/thread/threadpool.hpp>

class SThreadedNBList: public CSNeighbourList
{
public:
  SThreadedNBList(const XMLNode&, DYNAMO::SimData* const);

  SThreadedNBList(DYNAMO::SimData* const, CSSorter*, size_t threadCount);

  virtual void addEvents(const Particle&);
  
  virtual void operator<<(const XMLNode&);

  virtual void fullUpdate(const Particle& part);

  virtual void fullUpdate(const Particle& p1, const Particle& p2);

  void threadAddLocalEvent(const Particle& part, 
			   const size_t id,
			   magnet::thread::Mutex& sorterLock);

  void threadAddIntEvent(const Particle& part, 
			 const size_t id,
			 magnet::thread::Mutex& sorterLock);

  void spawnThreadAddLocalEvent1(const Particle& part, 
				const size_t& id);

  void spawnThreadAddLocalEvent2(const Particle& part, 
				 const size_t& id);

  void threadStreamParticles(const size_t id) const;

  void streamParticles(const Particle& part, const size_t& id) const;

  void addEvents2(const Particle& part, const size_t& id) const;

protected:
  virtual void outputXML(xml::XmlStream&) const;

  void addEventsInit(const Particle&);

  void addGlobal(const Particle& p1, const magnet::ClonePtr<Global>& glob, magnet::thread::Mutex& sorterLock);

  magnet::thread::ThreadPool _threadPool;

  magnet::thread::Mutex _P1SorterLock;
  magnet::thread::Mutex _P2SorterLock;
  
};
