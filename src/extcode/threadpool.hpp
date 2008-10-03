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
/*! \file threadpool.hpp
 * \brief Contains the definition of CThreadPool
 */

#ifndef THREADPOOL_H
#define THREADPOOL_H

#include <list>
#include <boost/thread/thread.hpp>
#include <boost/thread/condition.hpp>

//Just something to make calling a member function in a class using a member function pointer easy.
#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))

/*! \brief A class providing a pool of worker threads that will
 *   execute "tasks" pushed to it.
 *
 * Some helper classes are provided to call a member function of a
 * class as a task. These are functors.
 *
 * - CThreadPool::task_noarg - No arguments
 * - CThreadPool::task_1arg - 1 argument
 * - CThreadPool::task_2arg - 2 arguments
 * 
 * This class will also run in 0 thread mode, where the controlling
 * process will execute the tasks when it enters the CThreadPool::wait() function.
 *
 * This was adapted from an implementation by devguy.com
 * http://www.devguy.com/bb/viewtopic.php?p=1039
 *
 * \bug Check the licence on this class
 *
 */
class CThreadPool
{	
public:
  
  /*! \brief This class will turn a class and a member function
   * pointer into a functor.
   *
   * See the constructor for more details
   */
  template<class T>
  class task_noarg { 
  public:
    /*! An example of how to make a functor from this is     
     * task_noarg<CSimulation>(&A_Simulation_class, CSimulation::runSimulation)
     * 
     * \param nobj A reference to the class to be called
     * \param nfunc The member function pointer to call on the nobj class
     */
    task_noarg(T& nobj, void (T::*nfunc)()): 
      obj(&nobj), pfunc(nfunc) {}
        
    void operator()() 
    { CALL_MEMBER_FN(*obj,pfunc) (); }
    
  private:
    T *obj;
    void (T::*pfunc)();
  };

  /*! \brief This class will turn a class and a member function with 1
   * argument into a functor.
   */
  template<class T, class D>
  class task_1arg { 
  public:
    /*! An example of how to make a functor from this is     
     * task_noarg<CSimulation, size_t>(&A_Simulation_class, CSimulation::setSimID, 2)
     * 
     * \param nobj A reference to the class to be called
     * \param nfunc The member function pointer to call on the nobj class
     * \param D The single argument to the member function
     */
    task_1arg(T *nobj, void (T::*nfunc)(D), D data): 
      jobData(data), obj(nobj), pfunc(nfunc) {}

    void operator()() 
    { CALL_MEMBER_FN(*obj,pfunc) (jobData); }

  private:   
    D jobData;
    T *obj;
    void (T::*pfunc)(D);
  };

  /*! \brief This class will turn a class and a member function with 2
   * arguments into a functor.
   *
   * See task_1arg for more details
   */
  template<class T, class D1, class D2>
  class task_2arg { 
  public:
    task_2arg(T *nobj, void (T::*nfunc)(D1,D2),D1 data1, D2 data2): 
      jobData1(data1), jobData2(data2), obj(nobj), pfunc(nfunc) {}

    void operator()() 
    { CALL_MEMBER_FN(*obj,pfunc) (jobData1, jobData2); }

  private:   
    D1 jobData1;
    D2 jobData2;
    T *obj;
    void (T::*pfunc)(D1,D2);
  };

  /*! \brief A type that will allow generic holding of various (void) functors
   */
  typedef boost::function0<void> Functor_T;
  
  /*! \brief Default Constructor
   *
   * This initialises the pool to 0 threads
   */
  CThreadPool();
      
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
  { return m_threads.size(); }
    
  /*! \brief Set a task to be completed by the pool.
   *
   * \param threadfunc A functor which describes the task to complete.
   */
  inline void invoke(const Functor_T& threadfunc)
  {   
    boost::mutex::scoped_lock lock1(m_mutex);    

    m_needThread.notify_all();

    m_waitingFunctors.push_back(threadfunc);

    if (m_nextFunctor == m_waitingFunctors.end())
      --m_nextFunctor;

    lock1.unlock();
  }
  
  /*! \brief Destructor
   *
   * Join all threads in the pool and wait until they are terminated.
   */
  ~CThreadPool() throw();

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
  
  CThreadPool (const CThreadPool&);
  CThreadPool& operator = (const CThreadPool&);
  
  /*! \brief Current length of the functor/task queue. */
  size_t m_nQueueSize;

  /*! \brief List of functors/tasks left to be assigned to a thread. */
  std::list<Functor_T> m_waitingFunctors;

  /*! \brief Next functor in the list to be assigned. */
  std::list<Functor_T>::iterator m_nextFunctor;
  
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

  /*! \brief When this is true threads will terminate themselves.
   */
  bool m_bStop;

  //friend struct beginThreadFunc;
  
  /*! \brief Functor and entry point for a new thread so a thread can
   * access data in the owning CThreadPool.
   */
  struct beginThreadFunc
  {
    beginThreadFunc(CThreadPool& impl)
      : m_impl(impl)
    {}
    
    void operator() ()
    {
      m_impl.beginThread();
    }   
    
    /*! \brief A reference back to the original thread pool so we can
     * share mutexes in a object orientated program.
     */
    CThreadPool &m_impl;
  };
  
  /*! \brief Thread worker loop, called by the threads beginThreadFunc.
   */
  void beginThread() throw();

  /*! \brief Halt the threadpool and terminate all the threads.
   */
  void stop();
};
   

#endif
