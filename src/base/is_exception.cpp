/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "is_exception.hpp"
#include "is_stream_op.hpp"

namespace DYNAMO {
  
  exception::exception(int line, const char* file, const char* funcname) throw():
    message("")
  {
    message += "\nException thrown at [";
    message +=  boost::lexical_cast<std::string>(file);
    message += ":"; 
    message += boost::lexical_cast<std::string>(line); 
    message += "]";
#ifndef __DYNAMO_NO_EXCEPTION_FUNCTION
    message += "\nIn ";
    message += funcname;
#endif
    message += "\n";
  }
  
  const char *
  exception::what() const throw()
  {
    formattedMsg = searchReplace(message,"\n",boost::lexical_cast<std::string>("\n") 
				 + boost::lexical_cast<std::string>(IC_red) 
				 + boost::lexical_cast<std::string>("Exception") 
				 + boost::lexical_cast<std::string>(" :") 
				 + boost::lexical_cast<std::string>(IC_reset));

    return formattedMsg.c_str();
  }
}
