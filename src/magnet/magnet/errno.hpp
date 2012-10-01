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
#include <cerrno>
#include <cstring>
#include <magnet/exception.hpp>
#include <sstream>

namespace magnet {
  /*! \brief A reentrant, and C++ form of the strerror C function.
    
    \param errnum The error number to generate a string for.
   */
  inline std::string strerror(const int errnum)
  {
#if (_POSIX_C_SOURCE >= 200112L || _XOPEN_SOURCE >= 600) && ! _GNU_SOURCE
    std::vector<char> buf(128);

    int val = strerror_r(errnum, &buf[0], buf.size());

    while (val == ERANGE) 
      { 
	buf.resize(buf.size() * 2); 
	val = strerror_r(errnum, &buf[0], buf.size());
      }

    if (val == EINVAL) 
      {
	std::ostringstream os;
	os << "Unknown error number passed to strerror: " << errnum;
	return os.str();
      }

    if (val == -1) return std::string("Call to strerror_r failed");

    return std::string(buf.begin(), buf.end());
#else
    char buf[512];
    char* realbuf = strerror_r(errnum, buf, 512);
    return std::string(realbuf);
#endif
  }
}

