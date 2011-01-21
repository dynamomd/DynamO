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
/*! \file threadpool.hpp
 * \brief Contains the definition of ThreadPool
 */

#pragma once

#include <sstream>
#include <iostream>

#include <magnet/function/delegate.hpp>
#include <magnet/thread/taskQueue.hpp>
#include <magnet/memory/pool.hpp>

namespace magnet {
  namespace thread {
    /*! \brief A class providing a pool of worker threads that will
     *   execute "tasks" pushed to it.
     * 
     * This class will also run in 0 thread mode, where the controlling
     * process will execute the tasks when it enters the ThreadPool::wait() function.
     *
     */
    class ThreadPool : public TaskQueue
    {	
    private:
      bool _exception_flag;
      std::ostringstream _exception_data;
  
      ThreadPool (const ThreadPool&);
      ThreadPool& operator = (const ThreadPool&);
  
      /*! \brief This mutex is to control access to write that an exception has occurred.
       */
      magnet::thread::Mutex _exception_mutex;

      /*! \brief Triggered every time a thread becomes available, to
       * notify the mother thread stuck in the wait() function.
       */
      magnet::thread::Condition _threadAvailable_condition;

      /*! \brief Triggered to wake threads when jobs are added to the queue.
       */
      magnet::thread::Condition _need_thread_mutex;

      magnet::thread::ThreadGroup _threads;

      size_t _idlingThreads;
      size_t _threadCount;

      bool _stop_flag;

    public:  
      /*! \brief Default Constructor
       *
       * This initialises the pool to 0 threads
       */
      inline ThreadPool():
	_exception_flag(false),
	_idlingThreads(0),
	_stop_flag(false)
      {}
      
      /*! \brief Set the number of threads in the pool
       *
       * This creates the specified amount of threads to populate the
       * pool. When lowering the number of threads this pool kills ALL
       * threads, but waits for all their current tasks to complete first,
       * then repopulates the pool.
       */
      inline void setThreadCount(size_t x)
      { 
	if (x == _threadCount) return;
	
	if (x < _threadCount)
	  {
	    //Stop all threads as we're shrinking our thread pool size
	    stop();
	    //All threads are dead, reset the kill switch
	    _stop_flag = false;
	  }
	
	//Add the required number of threads
	for (size_t i=_threads.size(); i < x; ++i)
	  _threads.create_thread(function::Task::makeTask(&ThreadPool::beginThread, this));
	
	_threadCount = x;
      }

      /*! \brief The current number of threads in the pool */
      inline size_t getThreadCount() const { return _threadCount; }

      //Actual queuer
      inline void queueTask(function::Task* threadfunc)
      {
	TaskQueue::queueTask(threadfunc);
	_need_thread_mutex.notify_all();
      }
  
      /*! \brief Destructor
       *
       * Join all threads in the pool and wait until they are terminated.
       */
      inline ~ThreadPool() throw() { stop(); }

      /*! \brief Wait for all tasks to complete.
       *
       * If there are no threads in the pool then this function will
       * actually make the waiting/mother process perform the tasks.
       */
      inline void wait()
      {
	if (_threadCount)
	  {
	    //We are in threaded mode! Wait until all tasks are gone and all threads are idling
	    magnet::thread::ScopedLock lock1(_queue_mutex);      
	    while (!_waitingFunctors.empty() || (_idlingThreads != _threadCount))
	      _threadAvailable_condition.wait(lock1);
	  }
	else
	  {
	    //Non threaded mode
	    while (!_waitingFunctors.empty())
	      {
		(*_waitingFunctors.front())();
		delete _waitingFunctors.front();
		_waitingFunctors.pop();
	      }
	  }
	
	if (_exception_flag) 
	  M_throw() << "Thread Exception found while waiting for tasks/threads to finish"
		    << _exception_data.str();
      }

      inline const size_t& getIdleThreadCount() { return _idlingThreads; }
  
    private:
      /*! \brief Thread worker loop, called by the threads beginThreadFunc.
       */
      inline void beginThread()
      {
	try
	  {
	    magnet::thread::ScopedLock lock1(_queue_mutex);
	    
	    while (!_stop_flag)
	      {
		if (_waitingFunctors.empty())
		  {
		    ++_idlingThreads;
		    //Let whoever is waiting know that the thread is now available
		    _threadAvailable_condition.notify_all();
		    //And send it to sleep
		    _need_thread_mutex.wait(lock1);
		    --_idlingThreads;
		    continue;
		  }
		
		function::Task* func = _waitingFunctors.front();
		_waitingFunctors.pop();
		
		lock1.unlock();
		
		try { (*func)(); }
		catch(std::exception& cep)
		  {		  
		    //Mark the main process to throw an exception as soon as possible
		    magnet::thread::ScopedLock lock2(_exception_mutex);
		    
		    _exception_data << "\nTHREAD: Task threw an exception:-"
				     << cep.what();
		    
		    _exception_flag = true;
		  }
		
		delete func;
		
		lock1.lock();
	      }
	  }
	catch (std::exception& p)
	  {
	    std::cout << "\nTHREAD :Catastrophic Failure of thread!!! System will Hang, Aborting!"
		      << p.what();
	    throw;
	  }
      }
      
      /*! \brief Halt the threadpool and terminate all the threads.
       */
      inline void stop()
      {
	// _stop_flag must be set to true in a critical section.  Otherwise
	// it is possible for a thread to miss notify_all and never
	// terminate.
	{
	  magnet::thread::ScopedLock lock1(_queue_mutex);
	  _stop_flag = true;       
	}
	
	_need_thread_mutex.notify_all();
	_threads.join_all();
      }

    };
  }
}
