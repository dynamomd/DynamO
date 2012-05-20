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

#include <queue>
#include <magnet/thread/mutex.hpp>
#include <magnet/function/task.hpp>

namespace magnet {
  namespace thread {
    class TaskQueue
    {
    public:
      //Actual queuer
      inline virtual void queueTask(function::Task* threadfunc)
      {
	thread::ScopedLock lock1(_queue_mutex);    
	_waitingFunctors.push(threadfunc);
      }

      inline virtual void queueTasks(std::vector<function::Task*>& threadfuncs)
      {
	thread::ScopedLock lock1(_queue_mutex);

	for (std::vector<function::Task*>::const_iterator iPtr = threadfuncs.begin();
	     iPtr != threadfuncs.end(); ++iPtr)
	  _waitingFunctors.push(*iPtr);

	threadfuncs.clear();
      }

      void drainQueue()
      {
	_queue_mutex.lock();

	while (!_waitingFunctors.empty())
	  {
	    function::Task* task = _waitingFunctors.front();
	    _waitingFunctors.pop();
	    _queue_mutex.unlock();
	    (*task)();
	    delete task;
	    _queue_mutex.lock();
	  }
	_queue_mutex.unlock();
      }

      ~TaskQueue()
      {
	thread::ScopedLock lock1(_queue_mutex);

	//Just drain the queue
	while (!_waitingFunctors.empty())
	  {
	    delete _waitingFunctors.front();
	    _waitingFunctors.pop();
	  }
      }

    protected:
      std::queue<function::Task*> _waitingFunctors;
      thread::Mutex _queue_mutex;

    };
  }
}
