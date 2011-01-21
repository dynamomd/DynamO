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

#include <magnet/function/delegate.hpp>
//#include <magnet/memory/pool.hpp>

namespace magnet {
  namespace function {
    namespace detail { template<class T> struct typeWrapper { typedef T Type; }; }

    template<class T1> class Task0;
    template<class T1, class T2> class Task1;
    template<class T1, class T2, class T3> class Task2;
    template<class T1, class T2, class T3, class T4> class Task3;
    
    class Task //: public magnet::memory::PoolAllocated
    {
    public:
      virtual void operator()() = 0;
      virtual Task* Clone() = 0;

      template<typename retT>
      inline static Task* makeTask(retT (*funcPtr)()) 
      { 
	return new Task0<retT>(funcPtr); 
      }
      
      template<typename retT, typename classT>
      inline static Task* makeTask(retT (classT::*funcPtr)(), classT* classPtr)
      { 
	return new Task0<retT>(magnet::function::MakeDelegate(classPtr, funcPtr)); 
      }
      
      template<typename retT, typename classT>
      inline static Task* makeTask(retT (classT::*funcPtr)() const, classT* classPtr)
      { 
	return new Task0<retT>(magnet::function::MakeDelegate(classPtr, funcPtr)); 
      }
      
      //1 Argument
      template<typename retT, typename arg1T>
      inline static Task* makeTask(retT (*funcPtr)(arg1T), typename detail::typeWrapper<arg1T>::Type arg1)
      { return new Task1<retT,arg1T>(funcPtr, arg1); }
      
      template<typename retT, typename classT, typename arg1T>
      inline static Task* makeTask(retT (classT::*funcPtr)(arg1T), classT* classPtr, 
			    typename detail::typeWrapper<arg1T>::Type arg1)
      { 
	return new Task1<retT, arg1T>(magnet::function::MakeDelegate(classPtr, funcPtr), arg1); 
      }
      
      template<typename retT, typename classT, typename arg1T>
      inline static Task* makeTask(retT (classT::*funcPtr)(arg1T) const, classT* classPtr, 
			    typename detail::typeWrapper<arg1T>::Type arg1)
      { 
	return new Task1<retT, arg1T>(magnet::function::MakeDelegate(classPtr, funcPtr), arg1); 
      }
      
      //2 Argument
      template<typename retT, typename arg1T, typename arg2T>
      inline static Task*
      makeTask(retT (*funcPtr)(arg1T, arg2T), 
	       typename detail::typeWrapper<arg1T>::Type arg1, typename detail::typeWrapper<arg2T>::Type arg2)
      { 
	return new Task2<retT, arg1T, arg2T>(funcPtr, arg1, arg2); 
      }
      
      template<typename retT, typename classT, typename arg1T, typename arg2T>

      inline static Task*
      makeTask(retT (classT::*funcPtr)(arg1T, arg2T), classT* classPtr, 
	       typename detail::typeWrapper<arg1T>::Type arg1, typename detail::typeWrapper<arg2T>::Type arg2)
      { 
	return new Task2<retT, arg1T, arg2T>(magnet::function::MakeDelegate(classPtr, funcPtr), arg1, arg2); 
      }
      
      template<typename retT, typename classT, typename arg1T, typename arg2T>

      inline static Task*
      makeTask(retT (classT::*funcPtr)(arg1T, arg2T) const, classT* classPtr, 
	       typename detail::typeWrapper<arg1T>::Type arg1, typename detail::typeWrapper<arg2T>::Type arg2)
      { 
	return new Task2<retT, arg1T, arg2T>(magnet::function::MakeDelegate(classPtr, funcPtr), arg1, arg2); 
      }
      
      //3 Argument
      template<typename retT, typename arg1T, typename arg2T, typename arg3T>
      inline static Task*
      makeTask(retT (*funcPtr)(arg1T, arg2T, arg3T), 
	       typename detail::typeWrapper<arg1T>::Type arg1, 
	       typename detail::typeWrapper<arg2T>::Type arg2,
	       typename detail::typeWrapper<arg3T>::Type arg3)
      { 
	return new Task3<retT, arg1T, arg2T, arg3T>(funcPtr, arg1, arg2, arg3); 
      }
      
      template<typename retT, typename classT, typename arg1T, typename arg2T, typename arg3T>
      inline static Task*
      makeTask(retT (classT::*funcPtr)(arg1T, arg2T, arg3T), classT* classPtr, 
	       typename detail::typeWrapper<arg1T>::Type arg1, 
	       typename detail::typeWrapper<arg2T>::Type arg2,
	       typename detail::typeWrapper<arg3T>::Type arg3)
      { 
	return new Task3<retT, arg1T, arg2T, arg3T>(magnet::function::MakeDelegate(classPtr, funcPtr), 
						    arg1, arg2, arg3); 
      }
      
      template<typename retT, typename classT, typename arg1T, typename arg2T, typename arg3T>
      inline static Task*
      makeTask(retT (classT::*funcPtr)(arg1T, arg2T, arg3T) const, classT* classPtr, 
	       typename detail::typeWrapper<arg1T>::Type arg1, 
	       typename detail::typeWrapper<arg2T>::Type arg2,
	       typename detail::typeWrapper<arg3T>::Type arg3)
      { 
	return new Task3<retT, arg1T, arg2T, arg3T>(magnet::function::MakeDelegate(classPtr, funcPtr), 
						    arg1, arg2, arg3); 
      }

    };
    
    template<class T>
    class Task0 : public Task
    {
      typedef function::Delegate0<T> function;

    public:
      inline Task0(function delegate): _delegate(delegate) {}
      
      virtual void operator()() { _delegate(); }

      virtual Task* Clone() { return new Task0<T>(*this); }
      
    private:
      function _delegate;
    };
    
    template<class T, class T1>
    class Task1 : public Task
    {
      typedef function::Delegate1<T1, T> function;
    public:
      inline Task1(function delegate, T1 arg1): 
	_delegate(delegate), 
	_arg1(arg1) 
      {}
      
      virtual void operator()() { _delegate(_arg1); }

      virtual Task* Clone() { return new Task1<T, T1>(*this); }
      
    private:
      function _delegate;
      T1 _arg1;
    };
    
    template<class T, class T1, class T2>
    class Task2 : public Task
    {
      typedef function::Delegate2<T1, T2, T> function;
    public:
      inline Task2(function delegate, T1 arg1, T2 arg2): 
	_delegate(delegate), 
	_arg1(arg1),  
	_arg2(arg2) 
      {}
      
      virtual void operator()() { _delegate(_arg1, _arg2); }
      
      virtual Task* Clone() { return new Task2<T, T1, T2>(*this); }

    private:
      function _delegate;
      T1 _arg1;
      T2 _arg2;
    };
    
    template<class T, class T1, class T2, class T3>
    class Task3 : public Task
    {
      typedef function::Delegate3<T1, T2, T3, T> function;
    public:
      inline Task3(function delegate, T1 arg1, T2 arg2, T3 arg3): 
	_delegate(delegate), 
	_arg1(arg1),  
	_arg2(arg2),
	_arg3(arg3)
      {}
      
      virtual void operator()() { _delegate(_arg1, _arg2, _arg3); }

      virtual Task* Clone() { return new Task3<T, T1, T2, T3>(*this); }
      
    private:
      function _delegate;
      T1 _arg1;
      T2 _arg2;
      T3 _arg3;
    };
  }
}
