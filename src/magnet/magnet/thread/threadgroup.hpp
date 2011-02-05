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
#include <vector>

#include <magnet/exception.hpp>
#include <magnet/function/task.hpp>
#include <magnet/thread/thread.hpp>

namespace magnet {
  namespace thread {
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
