#include "PoolAllocator.hpp"

PoolManager PoolManager::pool;
boost::mutex PoolManager::poolLock;
