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
#include <magnet/string/searchreplace.hpp>

namespace magnet {
  namespace string {
    namespace detail {
      /*! \brief Class to track the indentation of a line.
       */
      class Indentor
      {
      public:
	Indentor(size_t factor = 2):_indentation(0), _factor(factor) {}

	Indentor& operator++() { ++_indentation; return *this; }
	Indentor& operator--() { --_indentation; return *this; }

	friend std::ostream& operator<<(std::ostream& os, const Indentor& id)
	{ 
	  for (size_t i(0); i < id._indentation * id._factor; ++i) 
	    os << ' '; 
	  return os;
	}

      protected:
	size_t _indentation;
	size_t _factor;
      };
    }
    /*! \brief Formats C code by adding whitespace.
     * \param in The string containing the C source code.
     * \returns A formatted version of the C code.
     */
    inline std::string 
    format_code(std::string in)
    { 
      std::ostringstream os;
      detail::Indentor indent;

      bool disable_linebreaks = false;

      for (std::string::const_iterator cPtr = in.begin();
	   cPtr != in.end(); ++cPtr)
	{
	  os << *cPtr;
	  switch (*cPtr) {
	  case '(': disable_linebreaks = true; break;
	  case ')': disable_linebreaks = false; break;
	  case ';': if (!disable_linebreaks) os << '\n' << indent; break;
	  case '{': os << '\n' << ++indent; break;
	  case '}': os << '\n' << --indent; break;
	  }
	}

      return os.str(); 
    }
  }
}
