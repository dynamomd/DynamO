/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

CThreadPool::CThreadPool():
  ExceptionThrown(false)
{
  m_bStop = false;
  m_nextFunctor = m_waitingFunctors.end();
}

void 
CThreadPool::setThreadCount(size_t x)
{ 
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
}

CThreadPool::~CThreadPool() throw()
{
  stop();
}

void 
CThreadPool::wait()
{
  if (m_threads.size())
    {
      //We are in threaded mode!
      boost::mutex::scoped_lock lock1(m_mutex);      
      while (!m_waitingFunctors.empty())
	{
	  m_threadAvailable.wait(lock1);
	}
      
      if (ExceptionThrown) 
	I_throw() << "Thread Exception found while waiting for tasks/threads to finish";
    }
  else
    {
      //Non threaded mode
      BOOST_FOREACH(Functor_T& afunctor, m_waitingFunctors)
	afunctor();
      
      m_waitingFunctors.clear();
    }
}

void 
CThreadPool::stop()
{
  // m_bStop must be set to true in a critical section.  Otherwise
  // it is possible for a thread to miss notify_all and never
  // terminate.
  boost::mutex::scoped_lock lock1(m_mutex);
  m_bStop = true;
  lock1.unlock();
  
  m_needThread.notify_all();
  m_threads.join_all();
}

void 
CThreadPool::beginThread() throw()
{
  try
    {
      boost::mutex::scoped_lock lock1(m_mutex);
      
      for(;;)
	{
	  if (m_bStop)
	    break;
	  
	  if (m_nextFunctor == m_waitingFunctors.end())
	    {
	      // Wait until someone needs a thread
	      m_needThread.wait(lock1);
	    }
	  else
	    {
	      std::list<Functor_T>::iterator iter = m_nextFunctor++;
	      
	      lock1.unlock();
	      
	      try
		{
		  (*iter)();
		}
	      catch(DYNAMO::Exception& cep)
		{
		  cep << "\nError in thread, task threw an exception, Handling gracefully";
		  std::cout << cep.what();
		  
		  //Mark the main process to throw an exception as soon as possible
		  boost::mutex::scoped_lock lock2(m_exception);
		  ExceptionThrown = true;
		}
	      
	      lock1.lock();
	      m_waitingFunctors.erase(iter);
	      lock1.unlock();
	      
	      m_threadAvailable.notify_all();
	      
	      lock1.lock();
	    }
	}
    }
  catch (DYNAMO::Exception& p)
    {
      p << "\nTHREAD :Catastrophic Failure of thread!!! System will Hang, Aborting!";
      std::cout << p.what();
      throw;
    }
  catch(...)
    {
      // This is a real problem.  Since this thread is about to die, m_nThreads
      // will be out of sync with the actual number of threads.  But this should
      // not occur except when something really bad happens.
      // not thread safe.. who cares
      I_throw() << "THREAD: Major Threading Error, unidentified exception";  
    }
}    
