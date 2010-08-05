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
#include "include/FastDelegate/FastDelegate.h"
#include "include/FastDelegate/FastDelegateBind.h"

#ifndef DYNAMO_CONDOR
#include <boost/thread/thread.hpp>
#include <boost/thread/condition.hpp>
#endif

#include "PoolAllocator.hpp"

/*! \brief A class providing a pool of worker threads that will
 *   execute "tasks" pushed to it.
 *
 * Some helper classes are provided to call a member function of a
 * class as a task. These are functors.
 *
 * - ThreadPool::task_noarg - No arguments
 * - ThreadPool::task_1arg - 1 argument
 * - ThreadPool::task_2arg - 2 arguments
 * 
 * This class will also run in 0 thread mode, where the controlling
 * process will execute the tasks when it enters the ThreadPool::wait() function.
 *
 * This was adapted from an implementation by devguy.com
 * http://www.devguy.com/bb/viewtopic.php?p=1039
 *
 * \bug Check the licence on this class
 * \bug Generalise the exception handling
 *
 */
class ThreadPool
{	
public:  
  class Task: public PoolAllocated
  {
  public:
    virtual void operator()() = 0;
  };

  template<class T>
  class Task0 : public Task
  {
    typedef fastdelegate::FastDelegate0<T> function;
  public:
    Task0(function delegate): _delegate(delegate) {}
    
    virtual void operator()() { _delegate(); }
    
  private:
    function _delegate;
  };

  template<class T, class T1>
  class Task1 : public Task
  {
    typedef fastdelegate::FastDelegate1<T1, T> function;
  public:
    Task1(function delegate, T1 arg1): 
      _delegate(delegate), 
      _arg1(arg1) 
    {}
    
    virtual void operator()() { _delegate(_arg1); }
    
  private:
    function _delegate;
    T1 _arg1;
  };

  template<class T, class T1, class T2>
  class Task2 : public Task
  {
    typedef fastdelegate::FastDelegate2<T1, T2, T> function;
  public:
    Task2(function delegate, T1 arg1, T2 arg2): 
      _delegate(delegate), 
      _arg1(arg1),  
      _arg2(arg2) 
   {}
    
    virtual void operator()() { _delegate(_arg1, _arg2); }
    
  private:
    function _delegate;
    T1 _arg1;
    T2 _arg2;
  };

  template<class T> struct typeWrapper { typedef T Type; };

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
#ifndef DYNAMO_CONDOR    
    return m_threads.size(); 
#else
    return 0;
#endif
  }
    
  /*! \brief Set a task to be completed by the pool.
   *
   * \param threadfunc A functor which describes the task to complete.
   */

  template<typename retT>
  inline void queue(retT (*funcPtr)()) 
  { queueTask(new Task0<retT>(funcPtr)); }

  template<typename retT, typename classT>
  void queue(retT (classT::*funcPtr)(), classT* classPtr)
  { queueTask(new Task0<retT>(fastdelegate::MakeDelegate(classPtr, funcPtr))); }

  //1 Argument
  template<typename retT, typename arg1T>
  inline void queue(retT (*funcPtr)(arg1T), typename typeWrapper<arg1T>::Type arg1)
  { queueTask(new Task1<retT,arg1T>(funcPtr, arg1)); }

  template<typename retT, typename classT, typename arg1T>
  inline void queue(retT (classT::*funcPtr)(arg1T), classT* classPtr, typename typeWrapper<arg1T>::Type arg1)
  { queueTask(new Task1<retT, arg1T>(fastdelegate::MakeDelegate(classPtr, funcPtr), arg1)); }

  template<typename retT, typename arg1T, typename arg2T>
  inline void queue(retT (*funcPtr)(arg1T, arg2T), 
		    typename typeWrapper<arg1T>::Type arg1, typename typeWrapper<arg2T>::Type arg2)
  { queueTask(new Task2<retT, arg1T, arg2T>(funcPtr, arg1, arg2)); }

  template<typename retT, typename classT, typename arg1T, typename arg2T>
  inline void queue(retT (classT::*funcPtr)(arg1T, arg2T), classT* classPtr, 
		    typename typeWrapper<arg1T>::Type arg1, typename typeWrapper<arg2T>::Type arg2)
  { queueTask(new Task2<retT, arg1T, arg2T>(fastdelegate::MakeDelegate(classPtr, funcPtr), arg1, arg2)); }

  inline void queueTask(Task* threadfunc)
  {
#ifndef DYNAMO_CONDOR   
    boost::mutex::scoped_lock lock1(m_mutex);    
    m_needThread.notify_all();
    m_waitingFunctors.push(threadfunc);

    lock1.unlock();
#else
    m_waitingFunctors.push(threadfunc);
#endif
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
  
private:
  /*! \brief Mark if an exception has thrown an exception so the
   * system can terminate gracefully.
   */
  bool ExceptionThrown;

  std::ostringstream ExceptionDetails;
  
  ThreadPool (const ThreadPool&);
  ThreadPool& operator = (const ThreadPool&);
  
  /*! \brief Current length of the functor/task queue. */
  size_t m_nQueueSize;

  /*! \brief List of functors/tasks left to be assigned to a thread. */
  std::queue<Task*> m_waitingFunctors;

#ifndef DYNAMO_CONDOR   
  /*! \brief A mutex to control access to job data.
   * 
   * This mutex is also used as part of the m_needThread condition and
   * also controls the m_bStop bool used to shutdown the queue. This
   * also protects m_waitingFunctors and m_nextFunctor.
   */
  boost::mutex m_mutex;

  /*! \brief This mutex is to control access to write that an exception has occurred.
   */
  boost::mutex m_exception;

  /*! \brief Triggered every time a thread becomes available, to
   * notify the mother thread stuck in the wait() function.
   */
  boost::condition m_threadAvailable;

  /*! \brief Triggered to wake threads when jobs are added to the queue.
   */
  boost::condition m_needThread;

  /*! \brief A collection of threads.
   */
  boost::thread_group m_threads;

  size_t _idlingThreads;

#endif

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
