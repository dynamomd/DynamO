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
#include <string>

namespace magnet {
  namespace string {
    /*! \brief Line wraps a string in a neat way.
     *
     * This function will wrap lines before they overflow
     *
     * \param in The string to be wrapped.
     * \param linelength The length of the line in characters.
     * \tparam hyphenate_long_words If true, any word longer than a line will be broken using hyphens
     * \returns The wrapped form of the in string.
     */
    template<bool hyphenate_long_words>
    inline std::string linewrap(std::string in, size_t linelength)
    {
      std::string::size_type linestart = -1, wordstart = -1;

      for (std::string::size_type curchar = 0; curchar <= in.size(); ++curchar)
	if ((in[curchar] == ' ')
	    || (curchar == in.size()))
	  {
	    //Will this word take the line length over its maximum?
	    if (curchar - linestart > linelength)
	      {
		//Check we can actually do a line break
		if (linestart != wordstart)
		  {
		    in[wordstart] = '\n';
		    linestart = wordstart;
		  }

		if (hyphenate_long_words)
		  while (curchar - linestart > linelength)
		    {
		      in.insert(linestart + linelength, std::string("\n-"));
		      linestart = linestart + linelength;
		      curchar += 2;
		    }
		else
		  if (curchar - linestart > linelength)
		    {
		      in[curchar] = '\n';
		      linestart = curchar;
		    }
	      }
	    wordstart = curchar;
	  }
	else if (in[curchar] == '\n')
	  wordstart = linestart = curchar;

      return in;
    }
  }
}
