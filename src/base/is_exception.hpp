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

#define I_throw() throw DYNAMO::Exception(__LINE__,__FILE__, __DYNAMO_EXCEPTION_FUNCTION)
#define I_throwRecoverable() throw DYNAMO::Exception(__LINE__,__FILE__, __DYNAMO_EXCEPTION_FUNCTION, true)

namespace DYNAMO
{
  class Exception
  {
  public:    
    inline ~Exception() throw() {}
    
    Exception(int line, const char* file, const char* funcname, bool recov = false) throw();
    
    template<class T>
    Exception& operator<<(T m) 
    { 
      message += boost::lexical_cast<std::string>(m);
      return *this;
    }
    
    std::string what() const throw();
    bool isRecoverable() const;
    
  private:
    std::string message;
    bool recoverable;
  }; 
}

#endif
