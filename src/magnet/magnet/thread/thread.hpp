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
#include <pthread.h>
#include <errno.h>

#include <magnet/exception.hpp>
#include <magnet/function/task.hpp>

namespace magnet {
  namespace thread {
    class Thread
    {
    public:
      Thread():_task(NULL), _joinable(false) {}
      
      ~Thread() 
      { if (_joinable) join(); }
      
      Thread(function::Task* task):
	_task(NULL), _joinable(false)
      { startTask(task); }

      inline void startTask(function::Task* task)
      {
	if (_joinable) join();
#ifdef MAGNET_DEBUG
        if (_task != NULL) M_throw() << "Task is not NULL!";
#endif
        _task = task;
        pthread_attr_t data;
        pthread_attr_init(&data);
        pthread_attr_setdetachstate(&data, PTHREAD_CREATE_JOINABLE);

        //We pass a pointer to the task
        if (pthread_create(&_thread, &data,
                           &threadEntryPoint,
                           static_cast<void*>(this)))
          M_throw() << "Failed to create a thread";
	
        pthread_attr_destroy(&data);
	_joinable = true;
      }
      
      inline void join()
      {
	int errval = pthread_join(_thread, NULL);
	_joinable = false;
	if (errval)
	  switch (errval)
	    {
	    case EINVAL:
	      M_throw() << "Failed to join thread, _task=" << _task
			<< "\n error is EINVAL";
	    case ESRCH:
	      M_throw() << "Failed to join thread, _task=" << _task
			<< "\n error is ESRCH";
	    case EDEADLK:
	      M_throw() << "Failed to join thread, _task=" << _task
			<< "\n error is EDEADLK";
	    default:
	      M_throw() << "Failed to join thread, _task=" << _task
			<< "\n error is UNKNOWN! errval = " << errval;
	    }
      }

      inline bool validTask() { return _task != NULL; }
      
    protected:
      //! Cannot be copied as the task pointer is a key communication variable to track thread activity
      Thread(const Thread&);
      //! See Thread(const Thread&)
      Thread& operator=(const Thread&);

      inline static void threadCleanup(void* arg) 
      {
	Thread* thread = reinterpret_cast<Thread*>(arg);
	delete thread->_task;
	thread->_task = NULL;
      }

      inline static void* threadEntryPoint(void* arg)
      {
        Thread* thread_ptr = reinterpret_cast<Thread*>(arg);
	pthread_cleanup_push(&Thread::threadCleanup, arg); //Set the cleanup function

	(*(*thread_ptr)._task)();

	pthread_cleanup_pop(1); //Pop the cleanup function and run it
	
	return NULL;
      }
      
      //! The pointer itself is volatile, but the object is not
      mutable function::Task* volatile _task;
      volatile bool _joinable;
      mutable pthread_t _thread;
    };
  }
}
