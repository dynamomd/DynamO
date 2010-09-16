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

#include <magnet/memory/pool.hpp>
#include <magnet/thread/mutex.hpp>

namespace magnet {
  namespace thread {
    class TaskQueue
    {
    public:
      //These tasks are allocated using a pool for speed
      class Task: public memory::PoolAllocated
      {
      public:
	virtual void operator()() = 0;
      };

      template<class T>
      class Task0 : public Task
      {
	typedef fastdelegate::FastDelegate0<T> function;
      public:
	inline Task0(function delegate): _delegate(delegate) {}
    
	virtual void operator()() { _delegate(); }
    
      private:
	function _delegate;
      };

      template<class T, class T1>
      class Task1 : public Task
      {
	typedef fastdelegate::FastDelegate1<T1, T> function;
      public:
	inline Task1(function delegate, T1 arg1): 
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
	inline Task2(function delegate, T1 arg1, T2 arg2): 
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

      template<class T, class T1, class T2, class T3>
      class Task3 : public Task
      {
	typedef fastdelegate::FastDelegate3<T1, T2, T3, T> function;
      public:
	inline Task3(function delegate, T1 arg1, T2 arg2, T3 arg3): 
	  _delegate(delegate), 
	  _arg1(arg1),  
	  _arg2(arg2),
	  _arg3(arg3)
	{}
    
	virtual void operator()() { _delegate(_arg1, _arg2, _arg3); }
    
      private:
	function _delegate;
	T1 _arg1;
	T2 _arg2;
	T3 _arg3;
      };

      template<class T> struct typeWrapper { typedef T Type; };

      template<typename retT>
      inline void queue(retT (*funcPtr)()) 
      { queueTask(new Task0<retT>(funcPtr)); }

      template<typename retT, typename classT>
      inline void queue(retT (classT::*funcPtr)(), classT* classPtr)
      { queueTask(new Task0<retT>(fastdelegate::MakeDelegate(classPtr, funcPtr))); }

      template<typename retT, typename classT>
      inline void queue(retT (classT::*funcPtr)() const, classT* classPtr)
      { queueTask(new Task0<retT>(fastdelegate::MakeDelegate(classPtr, funcPtr))); }

      //1 Argument
      template<typename retT, typename arg1T>
      inline void queue(retT (*funcPtr)(arg1T), typename typeWrapper<arg1T>::Type arg1)
      { queueTask(new Task1<retT,arg1T>(funcPtr, arg1)); }

      template<typename retT, typename classT, typename arg1T>
      inline void queue(retT (classT::*funcPtr)(arg1T), classT* classPtr, typename typeWrapper<arg1T>::Type arg1)
      { queueTask(new Task1<retT, arg1T>(fastdelegate::MakeDelegate(classPtr, funcPtr), arg1)); }

      template<typename retT, typename classT, typename arg1T>
      inline void queue(retT (classT::*funcPtr)(arg1T) const, classT* classPtr, typename typeWrapper<arg1T>::Type arg1)
      { queueTask(new Task1<retT, arg1T>(fastdelegate::MakeDelegate(classPtr, funcPtr), arg1)); }

      //2 Argument
      template<typename retT, typename arg1T, typename arg2T>
      inline void queue(retT (*funcPtr)(arg1T, arg2T), 
			typename typeWrapper<arg1T>::Type arg1, typename typeWrapper<arg2T>::Type arg2)
      { queueTask(new Task2<retT, arg1T, arg2T>(funcPtr, arg1, arg2)); }

      template<typename retT, typename classT, typename arg1T, typename arg2T>
      inline void queue(retT (classT::*funcPtr)(arg1T, arg2T), classT* classPtr, 
			typename typeWrapper<arg1T>::Type arg1, typename typeWrapper<arg2T>::Type arg2)
      { queueTask(new Task2<retT, arg1T, arg2T>(fastdelegate::MakeDelegate(classPtr, funcPtr), arg1, arg2)); }

      template<typename retT, typename classT, typename arg1T, typename arg2T>
      inline void queue(retT (classT::*funcPtr)(arg1T, arg2T) const, classT* classPtr, 
			typename typeWrapper<arg1T>::Type arg1, typename typeWrapper<arg2T>::Type arg2)
      { queueTask(new Task2<retT, arg1T, arg2T>(fastdelegate::MakeDelegate(classPtr, funcPtr), arg1, arg2)); }

      //3 Argument
      template<typename retT, typename arg1T, typename arg2T, typename arg3T>
      inline void queue(retT (*funcPtr)(arg1T, arg2T, arg3T), 
			typename typeWrapper<arg1T>::Type arg1, 
			typename typeWrapper<arg2T>::Type arg2,
			typename typeWrapper<arg3T>::Type arg3)
      { queueTask(new Task3<retT, arg1T, arg2T, arg3T>(funcPtr, arg1, arg2, arg3)); }

      template<typename retT, typename classT, typename arg1T, typename arg2T, typename arg3T>
      inline void queue(retT (classT::*funcPtr)(arg1T, arg2T, arg3T), classT* classPtr, 
			typename typeWrapper<arg1T>::Type arg1, 
			typename typeWrapper<arg2T>::Type arg2,
			typename typeWrapper<arg3T>::Type arg3)
      { queueTask(new Task3<retT, arg1T, arg2T, arg3T>(fastdelegate::MakeDelegate(classPtr, funcPtr), arg1, arg2, arg3)); }

      template<typename retT, typename classT, typename arg1T, typename arg2T, typename arg3T>
      inline void queue(retT (classT::*funcPtr)(arg1T, arg2T, arg3T) const, classT* classPtr, 
			typename typeWrapper<arg1T>::Type arg1, 
			typename typeWrapper<arg2T>::Type arg2,
			typename typeWrapper<arg3T>::Type arg3)
      { queueTask(new Task3<retT, arg1T, arg2T, arg3T>(fastdelegate::MakeDelegate(classPtr, funcPtr), arg1, arg2, arg3)); }


      void drainQueue()
      {
	thread::ScopedLock lock1(_mutex);

	while (!m_waitingFunctors.empty())
	  {
	    (*m_waitingFunctors.front())();
	    delete m_waitingFunctors.front();
	    m_waitingFunctors.pop();
	  }
      }

    protected:
      std::queue<Task*> _waitingFunctors;
      thread::Mutex _mutex;

      //Actual queuer
      inline void queueTask(Task* threadfunc)
      {
	thread::ScopedLock lock1(m_mutex);    
	m_waitingFunctors.push(threadfunc);
      }
    };
  }
}
