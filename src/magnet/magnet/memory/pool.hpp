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

#ifndef MAX_SMALL_OBJECT_SIZE
#define MAX_SMALL_OBJECT_SIZE 64
#endif

#include  <boost/pool/pool.hpp>
#include <magnet/thread/mutex.hpp>

namespace magnet {
  namespace memory {
    class PoolManager {
    public:
      /// singleton access

      inline static PoolManager& getPool()
      {
	static PoolManager pool;
	return pool;
      }

      inline static thread::Mutex& getLock()
      {  
	static thread::Mutex poolLock;
	return poolLock;
      }

      inline void* allocateMemory(size_t size) 
      {
	thread::ScopedLock lock(PoolManager::getLock());

	if (size > MAX_SMALL_OBJECT_SIZE)
	  return ::operator new(size);
    
	return m_pools[size - 1]->malloc();
      }

      inline void releaseMemory(void* deletable, size_t size) 
      {
	thread::ScopedLock lock(PoolManager::getLock());
	
	if (size > MAX_SMALL_OBJECT_SIZE) 
	  ::operator delete(deletable);
	else 
	  //use pool free. Don't delete null pointers
	  if (deletable)
	    m_pools[size - 1]->free(deletable);
      }

    private:

      inline PoolManager() {
	for (int i = 0; i < MAX_SMALL_OBJECT_SIZE; ++i)
	  m_pools[i] = new boost::pool<>(i + 1);
      }

      inline ~PoolManager() {
	for (int i = 0; i < MAX_SMALL_OBJECT_SIZE; ++i)
	  delete m_pools[i];
      }

      //no copying of the singleton
      PoolManager(const PoolManager&);
      const PoolManager& operator=(const PoolManager&);

      /// memory pool array. m_pools[n] corresponds to pool with objectSize==n+1.
      boost::pool<>* m_pools[MAX_SMALL_OBJECT_SIZE];

    };

    class PoolAllocated {
    public:
      inline static void* operator new(size_t size) {
	return PoolManager::getPool().allocateMemory(size);
      }

      inline static void operator delete(void* deletable, size_t size) {
	PoolManager::getPool().releaseMemory(deletable, size);
      }

      virtual ~PoolAllocated() {}
    };
  }
}
