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

//For the stack tracer
#ifdef MAGNET_DEBUG
# ifdef __GNUC__
#  include <execinfo.h>  // for backtrace
#  include <dlfcn.h>     // for dladdr
#  include <cxxabi.h>    // for __cxa_demangle
#  include <string>
#  include <sstream>
#  include <cstdio>
#  include <stdlib.h>
# endif
#endif

namespace magnet {

  /*! \brief This function will attempt to generate a string
   * representation of the stack at runtime.
   *
   * This functionality is only supported in the gcc compiler, if this
   * feature is unavailable an empty string is returned.
   *
   * \param skip The number of top stack levels to skip in the
   * trace. A skip of 1 will prevent the \ref stacktrace function from
   * appearing in the returned stack.
   */
#ifdef MAGNET_DEBUG
# ifdef __GNUC__
  inline std::string stacktrace(int skip = 1)
  {
    const int nMaxFrames = 128;
    void *callstack[nMaxFrames]; //The return addresses, including return offsets
    int nFrames = backtrace(callstack, nMaxFrames);
    char **symbols = backtrace_symbols(callstack, nFrames);
    
    std::ostringstream trace_buf;
    for (int i = skip; i < nFrames; i++) 
      {
	Dl_info info;
	if (dladdr(callstack[i], &info))
	  {
	    int status;
	    char *demangled = abi::__cxa_demangle(info.dli_sname, NULL, 0, &status);
	    
	    trace_buf << i << " " << callstack[i] << " "
		      << (status == 0 ? demangled : info.dli_sname)
		      << " + return offset=" << ((char *)callstack[i] - (char *)info.dli_saddr)
		      << "\n";
	    
	    free(demangled);
	  }
	else 
	  trace_buf << i << " " << callstack[i] << "\n";
      }
    free(symbols);

    if (nFrames == nMaxFrames)
      trace_buf << "[truncated]\n";

    return trace_buf.str();
  }
# else
  inline std::string stacktrace(int skip = 1) { return std::string(); }
# endif
#else
  inline std::string stacktrace(int skip = 1) { return std::string(); }
#endif
}
