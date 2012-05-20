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
#include <pthread.h>
#include <magnet/errno.hpp>
#include <magnet/function/task.hpp>

namespace magnet {
  namespace thread {
    class Mutex 
    {
    public:
      inline Mutex()
      {
#ifdef MAGNET_DEBUG
	pthread_mutexattr_t attr;
	int errornum = pthread_mutexattr_init(&attr);
	if (errornum)
	  M_throw() << "Could not initialize mutexattr\nError:" << magnet::strerror(errornum);
	
	errornum = pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_ERRORCHECK);
	if (errornum)
	  M_throw() << "Failed to set mutexattr type\nError:" << magnet::strerror(errornum);

	errornum = pthread_mutex_init(&_pthread_mutex, &attr);
	if (errornum)
	  M_throw() << "Failed to create the mutex\nError:" << magnet::strerror(errornum);

	errornum = pthread_mutexattr_destroy(&attr);
	if (errornum)
	  M_throw() << "Could not destroy mutexattr\nError:" << magnet::strerror(errornum);
#else
	int errornum = pthread_mutex_init(&_pthread_mutex, NULL);
	if (errornum)
	  M_throw() << "Failed to create the mutex\nError:" << magnet::strerror(errornum);
#endif
      }

      inline ~Mutex() 
      { 
	int errornum = pthread_mutex_destroy(&_pthread_mutex);
	if (errornum)
	  M_throw() << "Failed to destroy the mutex\nError:" << magnet::strerror(errornum);
      }

      inline void lock() 
      { 
	int errornum = pthread_mutex_lock(&_pthread_mutex);
	if (errornum)
	  M_throw() << "Failed to lock the mutex\nError:" << magnet::strerror(errornum);
     }

      inline bool try_lock() 
      { 
	int errornum = pthread_mutex_trylock(&_pthread_mutex);
	if ((errornum != 0) && (errornum != EBUSY))
	  M_throw() << "Failed to try_lock the mutex\nError:" << magnet::strerror(errornum);
	return !errornum;
      }

      inline void unlock() 
      {
	int errornum = pthread_mutex_unlock(&_pthread_mutex);
	if (errornum)
	  M_throw() << "Failed to unlock the mutex\nError:" << magnet::strerror(errornum);
      }

      inline pthread_mutex_t* native_handle()
      {
	return &_pthread_mutex;
      }

    protected:
      pthread_mutex_t _pthread_mutex;

    };

    class ScopedLock
    {
    public:
      inline ScopedLock(Mutex& mutex):
	_mutex(mutex),
	_is_locked(false)
      {
	lock();
      }

      inline ~ScopedLock() { if (is_locked()) _mutex.unlock(); }

      inline void unlock() { _mutex.unlock(); _is_locked = false; }

      inline void lock() { _mutex.lock(); _is_locked = true; }

      inline bool try_lock() { return _is_locked = _mutex.try_lock(); }

      inline operator const Mutex&() const { return _mutex; }
      inline operator Mutex&() { return _mutex; }

      inline bool is_locked() const { return _is_locked; }

    protected:
      Mutex& _mutex;
      bool _is_locked;
    };


    class Condition
    {
    public:
      inline Condition()
      {
	int errornum = pthread_cond_init(&_pthread_condition, NULL);
	if (errornum)
	  M_throw() << "Failed to initialise the condition variable\nError:" << magnet::strerror(errornum);
      }
      
      inline ~Condition()
      {
	int errornum = pthread_cond_destroy(&_pthread_condition);
	if (errornum)
	  M_throw() << "Failed to destroy the condition variable\nError:" << magnet::strerror(errornum);
      }

      inline void wait(Mutex& mutex)
      {
	int errornum = pthread_cond_wait(&_pthread_condition, mutex.native_handle());
	if (errornum)
	  M_throw() << "Failed to wait on the condition & mutex!\nError:" << magnet::strerror(errornum);
      }

      inline void notify_one()
      {
	int errornum = pthread_cond_signal(&_pthread_condition);
	if (errornum)
	  M_throw() << "Failed to notify_one() on the condition\nError:" << magnet::strerror(errornum);
      }

      inline void notify_all()
      {
	int errornum = pthread_cond_broadcast(&_pthread_condition);
	if (errornum)
	  M_throw() << "Failed to notify_all() on the condition\nError:" << magnet::strerror(errornum);
      }

    protected:
      pthread_cond_t _pthread_condition;
    };
  }
}
