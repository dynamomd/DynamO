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

#ifndef smrtPPointr_H
#define smrtPPointr_H

#include "../extcode/xmlwriter.hpp"
#include "../base/is_exception.hpp"

template<class T>
class smrtPlugPtr
{
 public:
  friend inline xmlw::XmlStream& operator<<(xmlw::XmlStream& XML, 
					    const smrtPlugPtr<T>& plugptr)
  { return XML << *(plugptr.obj); };
  

  inline void swap(smrtPlugPtr<T>& plug1)
  {
      std::swap(plug1.obj, obj);    
  }

  template<class A, class B>
  friend inline bool operator<(const smrtPlugPtr<A>&, 
			       const smrtPlugPtr<B>&);

  
  //Explicit stops the pointer being copied!
  //Assign a smart pointer
  inline explicit smrtPlugPtr(T* pointee) : obj(pointee) {}
  inline explicit smrtPlugPtr() : obj(NULL) {}

  //Copy the pointed to object as well
  inline smrtPlugPtr(const smrtPlugPtr<T> &p2):
    obj(NULL)
    {
      if (p2.obj != NULL)
	obj = p2->Clone();
    }

  inline smrtPlugPtr& operator=(const smrtPlugPtr<T>& p2)
    {
      if (obj != NULL)
	delete obj;

      obj = NULL;
      
      if (p2.obj != NULL)
	obj = p2->Clone();
	
      return *this;
    }

  //Delete the output plugin
  inline ~smrtPlugPtr() 
    { 
      if (obj != NULL) 
	delete obj; 
    }

  inline T* operator->()
    { 
#ifdef DYNAMO_DEBUG
      if (obj == NULL)
	I_throw() << "Attempting to dereference a null pluginpointer";
#endif
  
      return obj; 
    }
  
  inline T& operator*()
    { 
#ifdef DYNAMO_DEBUG
      if (obj == NULL)
	I_throw() << "Attempting to dereference a null pluginpointer";
#endif

      return *obj; 
    }

  inline T* release()
    { 
#ifdef DYNAMO_DEBUG
      if (obj == NULL)
	I_throw() << "Attempting to dereference a null pluginpointer";
#endif
      T* objptr = obj;
      obj = NULL;
      return objptr; 
    }

  inline const T* operator->() const
    { 
#ifdef DYNAMO_DEBUG
      if (obj == NULL)
	I_throw() << "Attempting to dereference a null pluginpointer";
#endif

      return obj; 
    }
  
  inline const T& operator*() const
    { 
#ifdef DYNAMO_DEBUG
      if (obj == NULL)
	I_throw() << "Attempting to dereference a null pluginpointer";
#endif
      
      return *obj; 
    }

  inline void set_ptr(T* A)
  {
    if (obj != NULL)
      delete obj;
    obj = A; 
  }
  
  inline T* get_ptr() { return obj; }

  inline const T* get_ptr() const { return obj; }
  
  inline bool operator!=(T* A) const { return (obj != A); }

  inline bool operator==(T* A) const { return (obj == A); }

 private:
   T* obj;
};

template<class A, class B>
inline bool operator<(const smrtPlugPtr<A>& p1, const smrtPlugPtr<B>& p2)
{ return (*p1.obj < *p2.obj); }

namespace std {
  template<class T>
  inline void swap(smrtPlugPtr<T>& plug1, smrtPlugPtr<T>& plug2)
  {
    plug1.swap(plug2);
  }
}

#endif
