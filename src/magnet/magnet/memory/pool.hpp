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

#ifndef MAX_SMALL_OBJECT_SIZE
#define MAX_SMALL_OBJECT_SIZE 64
#endif

#include <boost/pool/pool.hpp>
#include <magnet/thread/mutex.hpp>

namespace magnet {
  /*! \brief Namespace for memory management classes.*/
  namespace memory {    
    /*! \brief Namespace for memory management implementation
     * details.
     */
    namespace detail {
      /*! \brief Singleton class to manage the memory pools.
       *
       * This class manages several boost memory pools. Any classes
       * deriving from the \ref PoolAllocated class will use this
       * class to access memory pools to allocate their memory.
       * 
       */
      class PoolManager {
      public:
	/*! \brief Singleton access function. */
	inline static PoolManager& getPool()
	{
	  static PoolManager pool;
	  return pool;
	}
	
	/*! \brief Singleton pool lock access. */
	inline static thread::Mutex& getLock()
	{  
	  static thread::Mutex poolLock;
	  return poolLock;
	}
	
	/*! \brief Request some memory from a suitable pool. */
	inline void* allocateMemory(size_t size) 
	{
	  thread::ScopedLock lock(PoolManager::getLock());
	  
	  if (size > MAX_SMALL_OBJECT_SIZE)
	    return ::operator new(size);
	  
	  return m_pools[size - 1]->malloc();
	}
	
	/*! \brief Release some allocated memory from a suitable
	 * pool. 
	 */
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
	
	/*! \brief Hidden constructor as its a Singleton. */
	PoolManager(const PoolManager&);
	const PoolManager& operator=(const PoolManager&);
	
	/// memory pool array. m_pools[n] corresponds to pool with objectSize==n+1.
	boost::pool<>* m_pools[MAX_SMALL_OBJECT_SIZE];
      };
    }
    
    /*! \brief Base class for derived classes which want to be
     * allocated from a memory pool.
     *
     * Allocating objects from a memory pool is a way to speed up the
     * construction and deletion of small objects. The allocated
     * memory is not actually released but is instead cached for the
     * next allocation. This is handled by the \ref
     * detail::PoolManager.
     */
    class PoolAllocated {
    public:
      /*! \brief Specialized new operator to use the memory pool for
       * allocation. 
       */
      inline static void* operator new(size_t size) {
	return detail::PoolManager::getPool().allocateMemory(size);
      }

      /*! \brief Specialized delete operator to use the memory pool
       * for deallocation. 
       */
      inline static void operator delete(void* deletable, size_t size) {
	detail::PoolManager::getPool().releaseMemory(deletable, size);
      }
      
      virtual ~PoolAllocated() {}
    };
  }
}
