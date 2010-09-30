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

#include "threadpool.hpp"
#include <boost/foreach.hpp>
#include <magnet/exception.hpp>
#include <iostream>

namespace magnet {
  namespace thread {
    ThreadPool::ThreadPool():
      ExceptionThrown(false),
      _idlingThreads(0),
      m_bStop(false)
    {
      //  std::cerr << "\n m_mutex = " << m_mutex.native_handle()
      //	    << "\n m_exception = " << m_exception.native_handle()
      //	    << "\n m_threadAvailable = " << &m_threadAvailable  << " offset = " << (void*)((long)(&m_threadAvailable) + sizeof(pthread_mutex_t))
      //	    << "\n m_needThread = " << &m_needThread << " condition = " << (void*)((long)(&m_needThread) + sizeof(pthread_mutex_t))
      //	    << "\n";
    }

    void 
    ThreadPool::setThreadCount(size_t x)
    { 
      if (x == _threadCount) return;
  
      if (x < _threadCount)
	{
	  //Stop all threads as we're shrinking our thread pool size
	  stop();
	  //All threads are dead, reset the kill switch
	  m_bStop = false;
	}

      //Add the required number of threads
      for (size_t i=m_threads.size(); i < x; ++i)
	m_threads.create_thread(beginThreadFunc(*this));

      _threadCount = x;
    }

    ThreadPool::~ThreadPool() throw()
    {
      stop();
    }

    void 
    ThreadPool::wait()
    {
      if (_threadCount)
	{
	  //We are in threaded mode! Wait until all tasks are gone and all threads are idling
	  magnet::thread::ScopedLock lock1(_queue_mutex);      
	  while (!_waitingFunctors.empty() || (_idlingThreads != _threadCount))
	    m_threadAvailable.wait(lock1);
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

      {
	//This mutex is not needed here as the threads have stopped working!
	//boost::mutex::scoped_lock lock2(m_exception);

	if (ExceptionThrown) 
	  M_throw() << "Thread Exception found while waiting for tasks/threads to finish"
		    << ExceptionDetails.str();
      }
    }

    void 
    ThreadPool::stop()
    {
      // m_bStop must be set to true in a critical section.  Otherwise
      // it is possible for a thread to miss notify_all and never
      // terminate.
      {
	magnet::thread::ScopedLock lock1(_queue_mutex);
	m_bStop = true;       
      }

      m_needThread.notify_all();
      m_threads.join_all();
    }

    void 
    ThreadPool::beginThread() throw()
    {
      try
	{
	  magnet::thread::ScopedLock lock1(_queue_mutex);

	  while (!m_bStop)
	    {
	      if (_waitingFunctors.empty())
		{
		  ++_idlingThreads;
		  //Let whoever is waiting know that the thread is now available
		  m_threadAvailable.notify_all();
		  //And send it to sleep
		  m_needThread.wait(lock1);
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
		  magnet::thread::ScopedLock lock2(m_exception);
	      
		  ExceptionDetails << "\nTHREAD " << boost::this_thread::get_id() 
				   << ":Task threw an exception:-"
				   << cep.what();
	      
		  ExceptionThrown = true;
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
  }
}
