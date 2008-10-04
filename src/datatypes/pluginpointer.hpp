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
/*! \file pluginpointer.hpp
 *
 * \brief Holds the definition of the smrtPlugPtr class
 */


#ifndef smrtPPointr_H
#define smrtPPointr_H

#include "../extcode/xmlwriter.hpp"
#include "../base/is_exception.hpp"

/*! \brief A smart pointer with the ability to copy the polymorphic class it 
 * owns.
 *
 * This is just a simple smart pointer, to hold a heap allocated
 * object, and delete it when the class is destroyed. It has some
 * error checking to prevent dereferencing NULL pointers. More
 * importantly it has the ability to copy polymorphic classes when it
 * only contains a pointer to the base class.
 *
 * To copy a base class the class must contain a virtual function
 * Clone() returning heap allocated copied object.
 */
template<class T>
class smrtPlugPtr
{
 public:
  /*! \brief A helper function to allow contained objects to be
   * written to an xmlw::XmlStream.
   * 
   * \param XML xmlw::XmlStream to write the xml to.
   * \param plugptr smrtPlugPtr to write out to xml.
   * 
   * \return The xmlw::XmlStream so further writing can take place.
   */
  friend inline xmlw::XmlStream& operator<<(xmlw::XmlStream& XML, 
					    const smrtPlugPtr<T>& plugptr)
  { return XML << *(plugptr.obj); };
  

  /*! \brief A standard function to speed swapping pointers and pass
   * by val return statements.
   *
   * \param p1 The plugin with which to exhange the owned objects.
   */
  inline void swap(smrtPlugPtr<T>& p1)
  {
      std::swap(p1.obj, obj);    
  }

  /*! \brief A helper function to allow contained objects to be
   * compared.
   *   
   * This is used in comparing the CSystem events as these can update every 
   * event.
   * 
   * \param p1 The LHS of the < operator
   * \param p2 The RHS of the < operator
   *   
   * \return True if plug1 < plug2
   */
  template<class A, class B>
  friend inline bool operator<(const smrtPlugPtr<A>& p1, 
			       const smrtPlugPtr<B>& p2);

  
  /*! \brief Allows a smrtPlugPtr to be RAII
   *
   * \param pointee The newly allocated object to be owned by the smrtPlugPtr
   *
   * \warning The object passed as the pointee should be allocated
   * as:- smrtPlugPtr<BASETYPE> new_plug_ptr(new DERIVED_TYPE(....));
   * to prevent loose standard pointers.
   *
   */
  inline explicit smrtPlugPtr(T* pointee) : obj(pointee) {}

  /*! \brief The default constructor
   * 
   * Initialises obj to NULL which allows the error checking to work.
   */
  inline explicit smrtPlugPtr() : obj(NULL) {}

  /*! \brief The optimised copy constructor.
   * 
   * This copies the pointed to object using the virtual Clone()
   * method provided by the object. This should return a pointer to a
   * newly allocated object, which is then stored by its base pointer.
   *
   * This is optimal compared to the assignment operator as the class
   * doesn't have to check its empty first.
   *
   * \param p2 The smrtPlugPtr to be copied.
   */
  inline smrtPlugPtr(const smrtPlugPtr<T> &p2):
    obj(NULL)
    {
      if (p2.obj != NULL)
	obj = p2->Clone();
    }

  /*! \brief The assignment operator
   * 
   * This copies the pointed to object using the virtual Clone()
   * method provided by the object. This should return a pointer to a
   * newly allocated object, which is then stored by its base pointer.
   *
   * The plugin pointer must delete the current obj stored if it
   * exists.
   * 
   * \param p2 The smrtPlugPtr to be copied.
   */
  inline smrtPlugPtr& operator=(const smrtPlugPtr<T>& p2)
    {
      if (obj != NULL)
	delete obj;

      obj = NULL;
      
      if (p2.obj != NULL)
	obj = p2->Clone();
	
      return *this;
    }

  /*! \brief The Destructor.
   *
   * Deletes the obj stored if there is one.
   */
  inline ~smrtPlugPtr() 
    { 
      if (obj != NULL) 
	delete obj; 
    }
  
  /*! \brief Dereference the pointer, i believe this is recursive.
   */
  inline T* operator->()
    { 
#ifdef DYNAMO_DEBUG
      if (obj == NULL)
	D_throw() << "Attempting to dereference a null pluginpointer";
#endif
  
      return obj; 
    }
  
  /*! \brief Dereference the object.
   */
  inline T& operator*()
    { 
#ifdef DYNAMO_DEBUG
      if (obj == NULL)
	D_throw() << "Attempting to dereference a null pluginpointer";
#endif

      return *obj; 
    }

  /*! \brief Return a pointer the current obj stored and stop managing
      the obj.
   */
  inline T* release()
    { 
#ifdef DYNAMO_DEBUG
      if (obj == NULL)
	D_throw() << "Attempting to dereference a null pluginpointer";
#endif
      T* objptr = obj;
      obj = NULL;
      return objptr; 
    }

  /*! \brief Dereference the const pointer, i believe this is recursive.
   */
  inline const T* operator->() const
    { 
#ifdef DYNAMO_DEBUG
      if (obj == NULL)
	D_throw() << "Attempting to dereference a null pluginpointer";
#endif

      return obj; 
    }
  
  /*! \brief Dereference the const pointer, i believe this is recursive.
   */
  inline const T& operator*() const
    { 
#ifdef DYNAMO_DEBUG
      if (obj == NULL)
	D_throw() << "Attempting to dereference a null pluginpointer";
#endif
      
      return *obj; 
    }

  /*! \brief Change the stored obj.
   *
   * This function will delete the current obj stored if it exists.
   *
   * \param A New object to own.
   */
  inline void set_ptr(T* A)
  {
    if (obj != NULL)
      delete obj;
    obj = A; 
  }
  
  /*! \brief Returns the pointer currently stored! UNSAFE
   *
   * \bug Remove this function as its unsafe.
   *
   * \return Current object stored in obj
   */
  inline T* get_ptr() { return obj; }

  /*! \brief Returns the const pointer currently stored.
   *
   * \bug I believe this is memory safe as I don't think you can
   * delete a const *.
   *
   * \return Current object stored in obj.
   */
  inline const T* get_ptr() const { return obj; }
  
  
  /*! \brief Returns true if no object is stored by the smrtPlugPtr.
   */
  inline bool empty() const { return obj == NULL; }

 private:
  
  /*! \brief A pointer to the object owned by the smrtPlugPtr
   */
   T* obj;
};

/*! \brief Comparison operator to compare the contents of two smrtPlugPtr
 *
 * \param p1 LHS
 * \param p2 RHS
 *
 * \return True if contents of p1 < contents of p2
 */
template<class A, class B>
inline bool operator<(const smrtPlugPtr<A>& p1, const smrtPlugPtr<B>& p2)
{ return (*p1.obj < *p2.obj); }

namespace std {
  /*! \brief A friend function to speed swapping plugins and return by value.
   *
   * \param p1 Plugin to swap
   * \param p2 Plugin to swap
   *
   */
  template<class T>
  inline void swap(smrtPlugPtr<T>& p1, smrtPlugPtr<T>& p2)
  {
    p1.swap(p2);
  }
}

#endif
