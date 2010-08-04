#pragma once

#include <boost/pool/pool.hpp>
#include <boost/thread/thread.hpp>

#ifndef MAX_SMALL_OBJECT_SIZE
#define MAX_SMALL_OBJECT_SIZE 64
#endif

class PoolManager {
public:
  /// singleton access
  static PoolManager pool;
  static boost::mutex poolLock;

  void* allocateMemory(size_t size) {
    boost::mutex::scoped_lock lock(poolLock);

    if (size > MAX_SMALL_OBJECT_SIZE)
      return ::operator new(size);
    
    return m_pools[size - 1]->malloc();
  }

  void releaseMemory(void* deletable, size_t size) {
    boost::mutex::scoped_lock lock(poolLock);
    
    if (size > MAX_SMALL_OBJECT_SIZE) {
      ::operator delete(deletable);
    } else {
      //use pool free. Don't delete null pointers
      if (deletable)
	m_pools[size - 1]->free(deletable);
    }
  }

  ~PoolManager() {
    for (int i = 0; i < MAX_SMALL_OBJECT_SIZE; ++i)
      delete m_pools[i];
  }

  PoolManager() {
    for (int i = 0; i < MAX_SMALL_OBJECT_SIZE; ++i)
      m_pools[i] = new boost::pool<>(i + 1);
  }

private:

  //no copying of the singleton
  PoolManager(const PoolManager&);
  const PoolManager& operator=(const PoolManager&);

  /// memory pool array. m_pools[n] corresponds to pool with objectSize==n+1.
  boost::pool<>* m_pools[MAX_SMALL_OBJECT_SIZE];

};

class PoolAllocated {
public:
  static void* operator new(size_t size) {
    return PoolManager::pool.allocateMemory(size);
  }

  static void operator delete(void* deletable, size_t size) {
    PoolManager::pool.releaseMemory(deletable, size);
  }

  virtual ~PoolAllocated() {}
};
