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
      {
	if (other._obj != NULL)
	  {
	    other._mutex->lock();
	    _obj = other._obj;
	    ++(*(_counter = other._counter));
	    _mutex = other._mutex;	    
	    _mutex->unlock();
	  }
      }

      inline ~RefPtr() { release(); }
      
      inline bool isValid() const { return _obj != NULL; }

      inline RefPtr<T>& operator=(const RefPtr<T>& other)
      {
	release();

	if (other._obj != NULL)
	  {
	    other._mutex->lock();	    
	    _obj = other._obj;
	    _mutex = other._mutex;
	    _counter = other._counter;
	    ++(*_counter);
	    _mutex->unlock();
	  }

	return *this;
      }

      inline T& operator*() { return *_obj; }
      inline const T& operator*() const { return *_obj; }

      inline T* operator->() { return _obj; }
      inline const T* operator->() const { return _obj; }

      inline void release()
      {
	if (_obj != NULL)
	  {
	    _mutex->lock();
	    --(*_counter);
	    if(*_counter == 0)
	      {
		delete _obj;
		delete _counter;
		_mutex->unlock();
		delete _mutex;
	      }
	    else
	      _mutex->unlock();

	    _obj = NULL;
	    _counter = NULL;
	    _mutex = NULL;
	  }
      }

      template<class T2> 
      T2& as() { return dynamic_cast<T2&>(*_obj); }

      template<class T2> 
      T2& as() const { return dynamic_cast<const T2&>(*_obj); }

    private:
      T* _obj;
      size_t *_counter;
      Mutex *_mutex;
    };
  }
}
