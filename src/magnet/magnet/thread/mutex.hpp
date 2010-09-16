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
#include <pthread.h>
#include <magnet/exception.hpp>

namespace magnet {
  namespace thread {
    class Mutex 
    {
    public:

      inline Mutex() 
      { 
	if (pthread_mutex_init(&_pthread_mutex, NULL))
	  M_throw() << "Failed to lock the mutex";
      }

      inline ~Mutex() { pthread_mutex_destroy(&_pthread_mutex); }

      inline void lock() 
      { 
	if (pthread_mutex_lock(&_pthread_mutex)) 
	  M_throw() << "Failed to lock the mutex.";
     }

      inline void unlock() 
      { 
	if (pthread_mutex_unlock(&_pthread_mutex))
	  M_throw() << "Failed to unlock the mutex.";
      }

    protected:

      pthread_mutex_t _pthread_mutex;

    };

    class ScopedLock
    {
    public:
      inline ScopedLock(Mutex& mutex):
	_mutex(mutex)
      {
	_mutex.lock();
      }

      inline ~ScopedLock() { _mutex.unlock(); }

    protected:
      Mutex& _mutex;
    };
  }
}
