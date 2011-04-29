
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
#include <magnet/exception.hpp>

namespace xml { class XmlStream; }

namespace magnet {
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
  class ClonePtr
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
    friend inline ::xml::XmlStream& operator<<(::xml::XmlStream& XML, 
					       const ClonePtr<T>& plugptr)
    { return XML << *(plugptr._obj); };
  

    /*! \brief A standard function to speed swapping pointers and pass
     * by val return statements.
     *
     * \param p1 The plugin with which to exhange the owned objects.
     */
    inline void swap(ClonePtr<T>& p1)
    {
      std::swap(p1._obj, _obj);    
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
    friend inline bool operator<(const ClonePtr<A>& p1, 
				 const ClonePtr<B>& p2);

  
    /*! \brief Allows a smrtPlugPtr to be RAII
     *
     * \param pointee The newly allocated object to be owned by the smrtPlugPtr
     *
     * \warning The object passed as the pointee should be allocated
     * as:- smrtPlugPtr<BASETYPE> new_plug_ptr(new DERIVED_TYPE(....));
     * to prevent loose standard pointers.
     *
     */
    inline explicit ClonePtr(T* pointee) : _obj(pointee) {}

    /*! \brief The default constructor
     * 
     * Initialises obj to NULL which allows the error checking to work.
     */
    inline explicit ClonePtr() : _obj(NULL) {}

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
    inline ClonePtr(const ClonePtr<T> &p2):
      _obj(NULL)
    {
      if (p2._obj != NULL)
	_obj = p2->Clone();
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
    inline ClonePtr& operator=(const ClonePtr<T>& p2)
    {
      if (_obj != NULL)
	delete _obj;

      _obj = NULL;
      
      if (p2._obj != NULL)
	_obj = p2->Clone();
	
      return *this;
    }

    /*! \brief The Destructor.
     *
     * Deletes the obj stored if there is one.
     */
    inline ~ClonePtr() 
    { 
      if (_obj != NULL) 
	delete _obj; 
    }
  
    /*! \brief Dereference the pointer, i believe this is recursive.
     */
    inline T* operator->()
    { 
#ifdef dynamo_DEBUG
      if (_obj == NULL)
	M_throw() << "Attempting to dereference a null pluginpointer";
#endif
  
      return _obj; 
    }
  
    /*! \brief Dereference the object.
     */
    inline T& operator*()
    { 
#ifdef dynamo_DEBUG
      if (_obj == NULL)
	M_throw() << "Attempting to dereference a null pluginpointer";
#endif

      return *_obj; 
    }

    /*! \brief Return a pointer the current obj stored and stop managing
      the obj.
    */
    inline T* release()
    { 
#ifdef dynamo_DEBUG
      if (_obj == NULL)
	M_throw() << "Attempting to dereference a null pluginpointer";
#endif
      T* objptr = _obj;
      _obj = NULL;
      return objptr; 
    }

    /*! \brief Dereference the const pointer, i believe this is recursive.
     */
    inline const T* operator->() const
    { 
#ifdef dynamo_DEBUG
      if (_obj == NULL)
	M_throw() << "Attempting to dereference a null pluginpointer";
#endif

      return _obj; 
    }
  
    /*! \brief Dereference the const pointer, i believe this is recursive.
     */
    inline const T& operator*() const
    { 
#ifdef dynamo_DEBUG
      if (_obj == NULL)
	M_throw() << "Attempting to dereference a null pluginpointer";
#endif
      
      return *_obj; 
    }

    /*! \brief Change the stored obj.
     *
     * This function will delete the current obj stored if it exists.
     *
     * \param A New object to own.
     */
    inline void set_ptr(T* A)
    {
      if (_obj != NULL)
	delete _obj;
      _obj = A; 
    }
  
    /*! \brief Returns the pointer currently stored! UNSAFE
     *
     * \bug Remove this function as its unsafe.
     *
     * \return Current object stored in obj
     */
    inline T* get_ptr() { return _obj; }

    /*! \brief Returns the const pointer currently stored.
     *
     * \bug I believe this is memory safe as I don't think you can
     * delete a const *.
     *
     * \return Current object stored in obj.
     */
    inline const T* get_ptr() const { return _obj; }
  
  
    /*! \brief Returns true if no object is stored by the smrtPlugPtr.
     */
    inline bool empty() const { return _obj == NULL; }

  private:
  
    /*! \brief A pointer to the object owned by the smrtPlugPtr
     */
    T* _obj;
  };

  /*! \brief Comparison operator to compare the contents of two smrtPlugPtr
   *
   * \param p1 LHS
   * \param p2 RHS
   *
   * \return True if contents of p1 < contents of p2
   */
  template<class A, class B>
  inline bool operator<(const ClonePtr<A>& p1, const ClonePtr<B>& p2)
  { return (*p1._obj < *p2._obj); }  
}

namespace std {
  /*! \brief A friend function to speed swapping plugins and return by value.
   *
   * \param p1 Plugin to swap
   * \param p2 Plugin to swap
   *
   */
  template<class T>
  inline void swap(magnet::ClonePtr<T>& p1, magnet::ClonePtr<T>& p2)
  {
    p1.swap(p2);
  }
}
