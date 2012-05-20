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

#include <exception>
#include <sstream>
#include <magnet/stacktrace.hpp>

#define M_throw()							\
  throw magnet::exception(__LINE__,__FILE__, __func__, magnet::stacktrace())

namespace magnet
{
  /*! \brief An exception class that works like a stream object.
   
    This class is thrown using the M_throw() macro like so :-
    M_throw() << "My custom line error";
   
    Remember to catch the exception class by reference! This prevents
    type conversions to base classes and maintains the virtual
    overrides.

    To allow the combining of multiple exception types, all exceptions
    should inherit virtually from base exception types.
  */
  class exception : public virtual std::exception
  {
  public:
    
    inline ~exception() throw() {}
    
    /*! \brief Constructor called by the M_throw() macro.
     
      \param line The line number in the source file.
      \param file The name of the source file.
      \param funcname The name of the function throwing the exception.
    */
    exception(int line, const char* file, const char* funcname, std::string stackTrace)
      throw() :
      _stackTrace(stackTrace)
    {
      _message << "\nException thrown in " << funcname
	       << " " << file << ":" << line << "]\n";
    }

    /*! \brief Copy constructor.  

      Please catch exceptions by reference, they should not be copied
      except at the throw site.
    */
    exception(const exception&e):
      _stackTrace(e._stackTrace)
    { _message << e._message.str(); }


    /*! \brief The stream operator engine for the class
     */
    template<class T>
    exception& operator<<(const T& m) 
    { 
      _message << m;
      return *this;
    }
    
    /*! \brief Get the stored message from the class.
      
      This message is not reentrant and is annoyingly C like. In a
      threaded environment this is still thread safe provided the
      exception try/catch statements are local to the thread, which
      ThreadPool ensures.
    */
    const char* what() const throw()
    {
      _whatval = _message.str();
      
      if (!_stackTrace.empty())
	_whatval += std::string("\nStack trace follows\n")
	  + _stackTrace;
      
      return _whatval.c_str();
    }

  private:
    /*! \brief Stores the message of the exception.
     */
    std::ostringstream _message;
    std::string _stackTrace;
    mutable std::string _whatval;
  }; 
}

