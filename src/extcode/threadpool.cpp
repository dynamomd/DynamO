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
#include "../base/is_exception.hpp"
#include <iostream>

ThreadPool::ThreadPool():
  ExceptionThrown(false),
  _idlingThreads(0),
  m_bStop(false)
{}

void 
ThreadPool::setThreadCount(size_t x)
{ 
#ifndef DYNAMO_CONDOR   
  if (x == m_threads.size()) return;
  
  if (x > m_threads.size())
    {
      for (size_t i=m_threads.size(); i < x; ++i)
	m_threads.create_thread(beginThreadFunc(*this));

      return;
    }
    
  if (m_threads.size() != 0)
    {
      stop();
      //All threads are dead, reset the kill switch
      m_bStop = false;
    }      
  
  for (size_t i=0; i<x; i++)
    m_threads.create_thread(beginThreadFunc(*this));
#else
  D_throw() << "Cannot use the thread pool as threading is disabled";
#endif
}

ThreadPool::~ThreadPool() throw()
{
  stop();
}

void 
ThreadPool::wait()
{
#ifndef DYNAMO_CONDOR   
  if (m_threads.size())
    {
      //We are in threaded mode!
      {
	boost::mutex::scoped_lock lock1(m_mutex);      
	while (_idlingThreads != m_threads.size())
	  m_threadAvailable.wait(lock1);
      }
      

      {
	boost::mutex::scoped_lock lock2(m_exception);
	if (ExceptionThrown) 
	  D_throw() << "Thread Exception found while waiting for tasks/threads to finish"
		    << ExceptionDetails;
      }
    }
  else
#endif
    {
      std::cerr << "\nThreadpool IS OFFLINE";
      //Non threaded mode
      while (!m_waitingFunctors.empty())
	{
	  (*m_waitingFunctors.front())();
	  delete m_waitingFunctors.front();
	  m_waitingFunctors.pop();
	}

      {
	boost::mutex::scoped_lock lock2(m_exception);
	if (ExceptionThrown)
	  D_throw() << "Thread Exception found while waiting for tasks/threads to finish";
      }
    }
}

void 
ThreadPool::stop()
{
  // m_bStop must be set to true in a critical section.  Otherwise
  // it is possible for a thread to miss notify_all and never
  // terminate.
#ifndef DYNAMO_CONDOR   
  boost::mutex::scoped_lock lock1(m_mutex);
  m_bStop = true;
  lock1.unlock();
  
  m_needThread.notify_all();
  m_threads.join_all();
#else
  m_bStop = true;
#endif
}

void 
ThreadPool::beginThread() throw()
{
#ifndef DYNAMO_CONDOR   
  try
    {
      boost::mutex::scoped_lock lock1(m_mutex);

      for(;;)
	{
	  if (m_bStop) break;
	  
	  if (m_waitingFunctors.empty())
	    {
	      ++_idlingThreads;
	      //Let whoever is waiting know that the thread is now available
	      m_threadAvailable.notify_all();
	      //And send it to sleep
	      m_needThread.wait(lock1);
	      --_idlingThreads;
	      continue;
	    }

	  Task* func = m_waitingFunctors.front();
	  m_waitingFunctors.pop();
	  
	  lock1.unlock();
	  
	  try { (*func)(); }
	  catch(std::exception& cep)
	    {		  
	      //Mark the main process to throw an exception as soon as possible
	      boost::mutex::scoped_lock lock2(m_exception);
	      
	      ExceptionDetails << "\nTHREAD :Error in thread, task threw an exception:-"
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

#else
  D_throw() << "Cannot beginThread as threading is disabled";
#endif      
}    
