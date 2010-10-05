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
#include <errno.h>
#include <vector>

#include <magnet/exception.hpp>
#include <magnet/function/task.hpp>

namespace magnet {
  namespace thread {
    class Mutex 
    {
    public:
      inline Mutex()
      { 
	if (pthread_mutex_init(&_pthread_mutex, NULL))
	  M_throw() << "Failed to create the mutex";
      }

      inline ~Mutex() 
      { 
	if (pthread_mutex_destroy(&_pthread_mutex) != 0)
	  M_throw() << "Failed to destroy the mutex";
      }

      inline void lock() 
      { 
	if (pthread_mutex_lock(&_pthread_mutex))
	  M_throw() << "Failed to lock the mutex.";
     }

      inline bool try_lock() 
      { 
	int const val = pthread_mutex_trylock(&_pthread_mutex);
	if ((val != 0) && (val != EBUSY))
	  M_throw() << "Failed to try_lock the mutex.";
	return !val;
      }

      inline void unlock() 
      {
	if (pthread_mutex_unlock(&_pthread_mutex))
	  M_throw() << "Failed to unlock the mutex.";
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
	if (pthread_cond_init(&_pthread_condition, NULL))
	  M_throw() << "Failed to initialise the condition variable";
      }
      
      inline ~Condition()
      {
	if (pthread_cond_destroy(&_pthread_condition))
	  M_throw() << "Failed to destroy the condition variable";
      }

      inline void wait(Mutex& mutex)
      {
	if (pthread_cond_wait(&_pthread_condition,
			      mutex.native_handle()) != 0)
	  M_throw() << "Failed to wait on the condition & mutex!";
      }

      inline void notify_one()
      {
	if (pthread_cond_signal(&_pthread_condition) != 0)
	  M_throw() << "Failed to notify_one() on the condition";
      }

      inline void notify_all()
      {
	if (pthread_cond_broadcast(&_pthread_condition) != 0)
	  M_throw() << "Failed to notify_all() on the condition";
      }

    protected:
      pthread_cond_t _pthread_condition;
    };

    class Thread
    {
    public:
      Thread():_task(NULL) {}
      
      ~Thread() 
      {
	if (_task != NULL)
	  M_throw() << "Destroying an active thread!";
      }
      
      Thread(function::Task* task):
	_task(task)
      {
	pthread_attr_t data;
	pthread_attr_init(&data);
	pthread_attr_setdetachstate(&data, PTHREAD_CREATE_JOINABLE);

	if (pthread_create(&_thread, &data, 
			   &threadEntryPoint,
			   static_cast<void*>(_task)))
	  M_throw() << "Failed to create a thread";

	pthread_attr_destroy(&data);
      }
      
      inline void join()
      {
	if (_task == NULL)
	  M_throw() << "Cannot join, this thread had no task!";

	if (pthread_join(_thread, NULL))
	  M_throw() << "Failed to join thread, _task=" << _task;

	//Mark this thread as joined/ended
	_task = NULL;
      }

      Thread& operator=(const Thread& other)
      {
	if (_task != NULL)
	  M_throw() << "Trying to assign to a thread which is still valid!";

	_task = other._task;
	other._task = NULL;
	_thread = other._thread;

	return *this;
      }

      inline bool validTask() { return _task != NULL; }
      
    protected:      
      Thread(const Thread&);
      
      
      inline static void* threadEntryPoint(void* arg)
      {
	(*reinterpret_cast<function::Task*>(arg))();

	//Now delete the task
	delete reinterpret_cast<function::Task*>(arg);

	return NULL;
      }
      
      mutable function::Task* _task;
      
      mutable pthread_t _thread;
    };

    class ThreadGroup
    {
    public:
      inline ThreadGroup() {}

      inline ~ThreadGroup() 
      {
	for (std::vector<Thread*>::iterator iPtr = _threads.begin();
	     iPtr != _threads.end(); ++iPtr)
	  delete *iPtr;
      }

      inline void create_thread(function::Task* task)
      {
	_threads.push_back(new Thread(task));
      }

      inline size_t size() const { return _threads.size(); }
      
      inline void join_all()
      {
	for (std::vector<Thread*>::iterator iPtr = _threads.begin();
	     iPtr != _threads.end(); ++iPtr)
	  (*iPtr)->join();
      }

    protected:
      std::vector<Thread*> _threads;

      ThreadGroup(const ThreadGroup&);
      ThreadGroup& operator=(const ThreadGroup&);
    };
  }
}
