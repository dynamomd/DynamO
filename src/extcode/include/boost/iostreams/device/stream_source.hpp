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

#ifndef STREAM_SOURCE_HPP
#define STREAM_SOURCE_HPP

#include <boost/iostreams/concepts.hpp>
#include <algorithm>

namespace boost {
  namespace iostreams {

    template<class T>
    class stream_source : public source {

      T& underlying_stream;

    public:
      explicit stream_source(T& us):underlying_stream(us) {}
      
      std::streamsize read(char* s, std::streamsize n)
      {
	if (underlying_stream.eof())
	  return -1;

	underlying_stream.read(s,n);
	
	std::streamsize nread = underlying_stream.gcount();

	if (nread != 0)
	  {
	    return nread;
	  }
	else
	  {
	    return -1;
	  }
      }
      
    };
  }
}

#endif
