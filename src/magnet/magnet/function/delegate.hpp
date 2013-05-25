/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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
#include <functional>
#include <unordered_map>

namespace magnet 
{
  /* \brief A small and efficient delegate class.
     
     The purpose of this class is to homogenise all types of function
     pointers into a standard representation. The approach is the same
     as discussed in the link below:
     http://www.codeproject.com/Articles/136799/Lightweight-Generic-C-Callbacks-or-Yet-Another-Del
     The only difference here is that C++11 allows us to implement his
     approach without macros and in a very short space.
     
     It's not quite as fast as Don Clugston's delegates, but its as fast
     as possible without using compiler internals. The compiler does
     have the opportunity to optimise to Don's approach.
  */
  typedef std::pair<std::size_t, std::size_t> DelegateKey;

  template<typename... Args> class Delegate { Delegate() = delete; };
  template<typename RetType, typename... Args> class Delegate<RetType(Args...)>
  {
  public:
    template <RetType (*Fun_addr)(Args...)>
    static inline Delegate create()
    {
      Delegate d;
      d._this_ptr = 0;
      d._shunt_ptr = [](void*, Args... args){ return (*Fun_addr)(args...); };
      return d;
    }

    template <typename T, RetType (T::*Mem_fun_addr)(Args...)>
    static inline Delegate create(T* obj)
    {
      Delegate d;
      d._this_ptr = obj;
      d._shunt_ptr = [](void* obj, Args... args){ return (static_cast<T*>(obj)->*Mem_fun_addr)(args...); };
      return d;
    }

    template <typename T, RetType (T::*Mem_fun_addr)(Args...) const>
    static inline Delegate create(const T* obj)
    {
      Delegate d;
      d._this_ptr = const_cast<T*>(obj);
      d._shunt_ptr = [](void* obj, Args... args){ return (static_cast<const T*>(obj)->*Mem_fun_addr)(args...); };
      return d;
    }
  
    RetType operator()(Args... arguments) const { return (*_shunt_ptr)(_this_ptr, arguments...); }

    bool operator==(const Delegate& od) const { return (_this_ptr == od._this_ptr) && (_shunt_ptr == od._shunt_ptr); }
    bool operator!=(const Delegate& od) const { return !operator==(od); }

    operator DelegateKey() const
    {
      static_assert(sizeof(std::size_t)==sizeof(void*), "Error, object pointers are not the same size as size_t!");
      static_assert(sizeof(std::size_t)==sizeof(RetType (*)(void *, Args...)), "Error, global function pointers are not the same size as size_t!");
      return DelegateKey(reinterpret_cast<std::size_t>(_this_ptr), 
			 reinterpret_cast<std::size_t>(_shunt_ptr));
    }

  private:
    //The standard guarantees that a void pointer can store the address
    //of an object.
    void* _this_ptr;
    RetType (*_shunt_ptr)(void *, Args...);
  };

  namespace detail {
    struct Delegate_hash {
      std::size_t operator()(const magnet::DelegateKey& key) const {
	//Standard boost hash combining
	std::size_t hash1 = std::hash<std::size_t>()(key.first);
	std::size_t hash2 = std::hash<std::size_t>()(key.second);
	return hash1 +  0x9e3779b9 + (hash2<<6) + (hash2>>2);
      }
    };
  }

  struct Tracked {
    std::unordered_map<DelegateKey, Tracked*, detail::Delegate_hash> _tracked_connections;

    virtual void remove_tracked(DelegateKey key) 
    { _tracked_connections.erase(key); }
    
    virtual ~Tracked()
    {
      for (const auto& connection : _tracked_connections)
	connection.second->remove_tracked(connection.first);
    }
  };

  template <typename... Args> class Signal { Signal() = delete; };
  template <typename Ret_t, typename... Args> class Signal<Ret_t (Args...)>: Tracked
  {
    typedef Delegate<Ret_t(Args...)> Delegate_t;
    std::unordered_map<DelegateKey, Delegate_t, detail::Delegate_hash> _slots;

    template <typename T>
    void connect_sfinae(Delegate_t key, typename T::Tracked* instance)
    { 
      _tracked_connections[key] = instance;
      instance->_tracked_connections[key] = this;
    }

    template <typename T>
    void disconnect_sfinae(Delegate_t key, typename T::Tracked* instance)
    {
      _tracked_connections.erase(instance);
      instance->_tracked_connections.erase(this);
    }

    template <typename T> void connect_sfinae(...) {}
    template <typename T> void disconnect_sfinae(...) {}

    virtual void remove_tracked(DelegateKey key) {
      Tracked::remove_tracked(key);
      _slots.erase(key);
    }

  public:
    template <Ret_t (*Fun_addr)(Args... )>
    void connect()
    { 
      Delegate_t delegate = Delegate_t::template create<Fun_addr>();
      _slots[delegate] = delegate;
    }

    template <typename T, Ret_t (T::*Mem_fun_addr)(Args...)>
    void connect(T* obj)
    { 
      Delegate_t delegate = Delegate_t::template create<T, Mem_fun_addr>(obj);
      _slots[delegate] = delegate;
      connect_sfinae<T>(delegate, obj);
    }

    template <typename T, Ret_t (T::*Mem_fun_addr)(Args...) const>
    void connect(T* obj) 
    {
      Delegate_t delegate = Delegate_t::template create<T, Mem_fun_addr>(obj);
      _slots[delegate] = delegate;
      connect_sfinae<T>(delegate, obj);
    }

    template <Ret_t (*Fun_addr)(Args... )>
    void disconnect()
    { 
      Delegate_t delegate = Delegate_t::template create<Fun_addr>();
      _slots.erase(delegate);
    }

    template <typename T, Ret_t (T::*Mem_fun_addr)(Args...)>
    void disconnect(T* obj)
    { 
      Delegate_t delegate = Delegate_t::template create<T, Mem_fun_addr>(obj);
      _slots.erase(delegate);
      disconnect_sfinae<T>(delegate, obj);
    }

    template <typename T, Ret_t (T::*Mem_fun_addr)(Args...) const>
    void disconnect(T* obj) 
    {
      Delegate_t delegate = Delegate_t::template create<T, Mem_fun_addr>(obj);
      _slots.erase(delegate);
      disconnect_sfinae<T>(delegate, obj);
    }

    void operator() (Args... args) const
    { for (const auto& slot : _slots) slot.second(args...); }
  };
}
