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

#ifndef THREADPOOL_H
#define THREADPOOL_H

// See http://www.devguy.com/bb/viewtopic.php?p=1039

#include <list>
#include <boost/thread/thread.hpp>
#include <boost/thread/condition.hpp>

#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))

template<class T>
class task_noarg { 
 public:
  task_noarg(T *nobj, void (T::*nfunc)()): 
    obj(nobj), pfunc(nfunc) {}

  void operator()() 
   {
     CALL_MEMBER_FN(*obj,pfunc) ();
   }

 private:   
  T *obj;
  void (T::*pfunc)();
};

template<class T, class D>
class task_1arg { 
 public:
  task_1arg(T *nobj, void (T::*nfunc)(D), D data): 
    jobData(data), obj(nobj), pfunc(nfunc) {}

  void operator()() 
   {
     CALL_MEMBER_FN(*obj,pfunc) (jobData);
   }

 private:   
  D jobData;
  T *obj;
  void (T::*pfunc)(D);
};

template<class T, class D1, class D2>
class task_2arg { 
 public:
  task_2arg(T *nobj, void (T::*nfunc)(D1,D2),D1 data1, D2 data2): 
    jobData1(data1), jobData2(data2), obj(nobj), pfunc(nfunc) {}

  void operator()() 
   {
     CALL_MEMBER_FN(*obj,pfunc) (jobData1, jobData2);
   }

 private:   
  D1 jobData1;
  D2 jobData2;
  T *obj;
  void (T::*pfunc)(D1,D2);
};

class CThreadPool
{	
 public:  
  typedef boost::function0<void> Functor_T;
  
  CThreadPool();
  
  explicit CThreadPool(size_t);
    
  void setMaxThreads(size_t);
  
  size_t getMaxThreads() const
  {
    return m_nMaxThreads;
  }
  
  // Register a function to be called each time a new thread is created, and run by that thread
  void setThreadCreatedListener(const Functor_T &f)
  {
    m_threadCreated = f;
  }
  
  // Call a function in a separate thread managed by the pool
  inline void invoke(const Functor_T& threadfunc)
  {   
    boost::mutex::scoped_lock lock1(m_mutex);    
    addFunctor(threadfunc);
    lock1.unlock();
    m_needThread.notify_all();
  }
  
  ~CThreadPool() throw();

  // This method is not thread-safe
  void stop();
  
  // This method is not thread-safe, only call by the mother thread
  void wait();
  
 private:
  bool ExceptionThrown;
  bool recoverable;
  
  CThreadPool (const CThreadPool&);
  CThreadPool& operator = (const CThreadPool&);
  
  size_t m_nQueueSize,
    m_nMaxThreads;  

  typedef std::list<Functor_T> Container_T;
  Container_T m_waitingFunctors;
  Container_T::iterator m_nextFunctor;
  Functor_T	m_threadCreated;
  
  boost::mutex m_mutex;
  boost::mutex m_exception;

  boost::condition 	m_threadAvailable, // triggered when a thread is available
    m_needThread;      // triggered when a thread is needed
  boost::thread_group m_threads;
  bool m_bStop;
  
  friend struct beginThreadFunc;
  
  void addFunctor(const Functor_T&);
  
  struct NoOp
  {
    void operator () () const
    {}
  };
  
  struct beginThreadFunc
  {
    beginThreadFunc(CThreadPool& impl)
      : m_impl(impl)
    {}
    
    void operator() ()
    {
      m_impl.beginThread();
    }   
    
    CThreadPool &m_impl;
  };
  
  // Thread entry point.  This method runs once per thread.
  void beginThread() throw();
};
   

#endif
