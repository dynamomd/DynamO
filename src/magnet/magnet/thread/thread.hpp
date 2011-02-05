/*  DYNAMO:- Event driven molecular dynamics simulator 
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
  }
}
