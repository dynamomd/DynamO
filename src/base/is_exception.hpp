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

#ifndef IS_Exception_H
#define IS_Exception_H

#include <exception>
#include <features.h>
#include <boost/lexical_cast.hpp>

# ifndef DOXYGEN_SHOULD_IGNORE_THIS
#  if defined __cplusplus ? __GNUC_PREREQ (2, 6) : __GNUC_PREREQ (2, 4)
#    define __DYNAMO_EXCEPTION_FUNCTION	__PRETTY_FUNCTION__
#  else
#   if defined __STDC_VERSION__ && __STDC_VERSION__ >= 199901L
#    define __DYNAMO_EXCEPTION_FUNCTION	__func__
#   else
#    define __DYNAMO_EXCEPTION_FUNCTION 'Unknown Function'
#    define __DYNAMO_NO_EXCEPTION_FUNCTION 1
#   endif
#  endif
# endif

#define D_throw() throw DYNAMO::exception(__LINE__,__FILE__, __DYNAMO_EXCEPTION_FUNCTION)

namespace DYNAMO
{
  /*! \brief An exception class that works like a stream object.
   *
   * This class is thrown using the D_throw() macro like so :-
   * D_throw() << "My custom line error";
   *
   * Remember to catch the exception class by reference!
   */
  class exception : public std::exception
  {
  public:
    
    inline ~exception() throw() {}
    
    /*! \brief Constructor called by the D_throw() macro.
     *
     * \param line The line number in the source file.
     * \param file The name of the source file.
     * \param funcname The name of the function throwing the exception.
     */
    exception(int line, const char* file, const char* funcname) throw();
    
    /*! \brief The stream operator engine for the class
     */
    template<class T>
    exception& operator<<(const T& m) 
    { 
      message += boost::lexical_cast<std::string>(m);
      return *this;
    }
    
    /*! \brief Get the stored message from the class.
     * 
     * This message is not reentrant and is annoyingly C like. In a
     * threaded environment this is still thread safe provided the
     * exception try/catch statements are local to the thread, which
     * ThreadPool ensures.
     * 
     * \bug Fix the non reentrant behaviour.
     */
    const char* what() const throw();

  private:
    /*! \brief Stores the message of the exception.
     */
    std::string message;

    /*! \brief Stores the formatted message returned by what()
     */
    mutable std::string formattedMsg;
  }; 
}

#endif
