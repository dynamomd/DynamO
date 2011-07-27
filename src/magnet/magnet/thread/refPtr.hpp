/*  dynamo:- Event driven molecular dynamics simulator 
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

#include <magnet/thread/mutex.hpp>

namespace magnet {
  namespace thread {
    template<class T>
    class RefPtr
    {
    public:
      inline RefPtr():
	_obj(NULL),
	_counter(NULL),
	_mutex(NULL)
      {}

      inline RefPtr(T* obj):
	_obj(obj),
	_counter(new size_t(1)),
	_mutex(new Mutex)
      {}

      inline RefPtr(const RefPtr<T>& other):
	_obj(NULL),
	_counter(NULL),
	_mutex(NULL)	
      { operator=(other); }

      inline ~RefPtr() { release(); }
      
      inline bool isValid() const { return _obj != NULL; }

      //template<class T2>
      inline RefPtr<T>& operator=(const RefPtr<T>& other)
      {
	release();
	if (other._obj == NULL) return *this;
	other._mutex->lock();
	_obj = other._obj;
	_mutex = other._mutex;
	_counter = other._counter;
	++(*_counter);
	_mutex->unlock();

	return *this;
      }

      inline T& operator*() { checkValid(); return *_obj; }
      inline const T& operator*() const { checkValid(); return *_obj; }

      inline T* operator->() { checkValid(); return _obj; }
      inline const T* operator->() const { checkValid(); return _obj; }

      inline void release()
      {
	if (_obj == NULL) return;

	_mutex->lock();
	--(*_counter);
	if(*_counter == 0)
	  {
	    delete _obj;
	    delete _counter;
	    _obj = NULL;
	    _counter = NULL;

	    _mutex->unlock();
	    delete _mutex;
	    _mutex = NULL;
	  }
	else
	  _mutex->unlock();	
      }

      template<class T2> 
      T2& as() { checkValid(); return dynamic_cast<T2&>(*_obj); }

      template<class T2> 
      T2& as() const { checkValid(); return dynamic_cast<const T2&>(*_obj); }

      template<class T2>
      inline bool operator==(const T2& other) const
      { return *_obj == other; }

      template<class T2>
      inline bool operator==(const RefPtr<T2>& other) const
      { return *_obj == *other; }

    private:
      void checkValid() const
      {
#ifdef MAGNET_DEBUG
	if (_obj == NULL) 
	  M_throw() << "Bad operation on invalid RefPtr";
#endif
      }
      
      T* _obj;
      size_t *_counter;
      Mutex *_mutex;
    };
  }
}
