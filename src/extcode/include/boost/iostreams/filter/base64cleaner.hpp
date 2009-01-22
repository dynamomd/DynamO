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
#ifndef BASE64CLEANER_HPP
#define BASE64CLEANER_HPP

#include <boost/iostreams/char_traits.hpp> // EOF, WOULD_BLOCK
#include <boost/iostreams/concepts.hpp>    // input_filter
#include <boost/iostreams/operations.hpp>  // get

namespace boost { namespace iostreams {

    class base64cleaner_input_filter : public input_filter {
    public:
      explicit base64cleaner_input_filter()
        : endofstream(false)
      { }

      template<typename Source>
      int get(Source& src)
      {
	if (endofstream) return EOF;

        int c;
        while (true) {
	  if ((c = boost::iostreams::get(src)) == EOF || c == WOULD_BLOCK)
	    break;

	  if (c =='<')
	    {
	      endofstream=true;
	      return EOF;
	    }

	  if ((c != '\n') && (c != ' '))
	    break;
        }
        return c;
      }

      template<typename Source>
      void close(Source&) { endofstream = false; }
    private:
      bool endofstream;
    };

  } 
} 

#endif
