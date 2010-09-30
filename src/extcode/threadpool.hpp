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
/*! \file threadpool.hpp
 * \brief Contains the definition of ThreadPool
 */

#pragma once

#include <queue>
#include <sstream>

#include <boost/thread/thread.hpp>

#include <magnet/function/delegate.hpp>
#include <magnet/thread/taskQueue.hpp>
#include <magnet/memory/pool.hpp>

/*! \brief A class providing a pool of worker threads that will
 *   execute "tasks" pushed to it.
 * 
 * This class will also run in 0 thread mode, where the controlling
 * process will execute the tasks when it enters the ThreadPool::wait() function.
 *
 * This was adapted from an implementation by devguy.com
 * http://www.devguy.com/bb/viewtopic.php?p=1039
 *
 */
namespace magnet {
  namespace thread {
    class ThreadPool : public TaskQueue
    {	
    public:  
      /*! \brief Default Constructor
       *
       * This initialises the pool to 0 threads
       */
      ThreadPool();
      
      /*! \brief Set the number of threads in the pool
       *
       * This creates the specified amount of threads to populate the
       * pool. When lowering the number of threads this pool kills ALL
       * threads, but waits for all their current tasks to complete first,
       * then repopulates the pool.
       */
      void setThreadCount(size_t);
  
      /*! \brief The current number of threads in the pool */
      size_t getThreadCount() const
      {
	return _threadCount; 
      }

      //Actual queuer
      inline void queueTask(function::Task* threadfunc)
      {
	TaskQueue::queueTask(threadfunc);
	m_needThread.notify_all();
      }
  
      /*! \brief Destructor
       *
       * Join all threads in the pool and wait until they are terminated.
       */
      ~ThreadPool() throw();

      /*! \brief Wait for all tasks to complete.
       *
       * If there are no threads in the pool then this function will
       * actually make the waiting/mother process perform the tasks.
       */
      void wait();

      const size_t& getIdleThreadCount() { return _idlingThreads; }
  
    private:
      /*! \brief Mark if an exception has thrown an exception so the
       * system can terminate gracefully.
       */
      bool ExceptionThrown;

      std::ostringstream ExceptionDetails;
  
      ThreadPool (const ThreadPool&);
      ThreadPool& operator = (const ThreadPool&);
  
      /*! \brief This mutex is to control access to write that an exception has occurred.
       */
      magnet::thread::Mutex m_exception;

      /*! \brief Triggered every time a thread becomes available, to
       * notify the mother thread stuck in the wait() function.
       */
      magnet::thread::Condition m_threadAvailable;

      /*! \brief Triggered to wake threads when jobs are added to the queue.
       */
      magnet::thread::Condition m_needThread;

      /*! \brief A collection of threads.
       */
      boost::thread_group m_threads;

      size_t _idlingThreads;
      size_t _threadCount;

      /*! \brief When this is true threads will terminate themselves.
       */
      bool m_bStop;

      //friend struct beginThreadFunc;
  
      /*! \brief Functor and entry point for a new thread so a thread can
       * access data in the owning ThreadPool.
       */
      struct beginThreadFunc
      {
	beginThreadFunc(ThreadPool& impl)
	  : m_impl(impl)
	{}
    
	void operator() ()
	{
	  m_impl.beginThread();
	}   
    
	/*! \brief A reference back to the original thread pool so we can
	 * share mutexes in a object orientated program.
	 */
	ThreadPool &m_impl;
      };
  
      /*! \brief Thread worker loop, called by the threads beginThreadFunc.
       */
      void beginThread() throw();

      /*! \brief Halt the threadpool and terminate all the threads.
       */
      void stop();
    };
  }
}
