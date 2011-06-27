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
#include <magnet/string/linewrap.hpp>
#include <iostream>
#include <sstream>

namespace magnet
{
  namespace stream {
    /*! \brief This class provides an std::ostream which wraps another
     * ostream but adds automatic formatting.
     *
     * The primary purpose of this class is to provide formatted
     * output for classes. The idea is that the output of each class
     * might be prefixed with some identifying information, and long
     * lines will be automatically wrapped. To acheive this, the end
     * of every chunk of information passed to this class must finish
     * with a std::endl. e.g.
     *
     * \code stream::FormattedOStream os(...);
     * os << "Some long text as part of a block of output, plus a number " << 20 
     *    << "but always finished with a endl." << std::endl; \endcode
     */
    class FormattedOStream : public std::ostream
    {
      
      /*! \brief A stream buffer which formats strings when flushed.
       *
       * This stream buffer overrides std::stringbuf so that when a flush is
       * called on the stream, the string buffer is formatted before output.
       */
      class FormatingStreamBuf: public std::stringbuf
      {

	/*! \brief The final destination of the formatted output..*/
        std::ostream&  _output;
	
	/*! \brief The maximum length of a formatted line before it is
	 *   wrapped. */
	size_t _linelength;
	
	/*! \brief Name to insert after newlines.*/
	std::string _prefix;

      public:
	FormatingStreamBuf(const std::string & prefix,
			   std::ostream& ostream,
			   const size_t linelength):
	  _output(ostream),
	  _linelength(linelength - prefix.size()),
	  _prefix(prefix)
	{}
	
	/*! \brief Copy constructor. */
	FormatingStreamBuf(const FormatingStreamBuf& os):
	  _output(os._output),
	  _linelength(os._linelength),
	  _prefix(os._prefix)
	{}

	std::ostream& getOutputStream() { return _output; }

	/*! \brief sync function override which actually performs the
	 *  output formatting before a std::flush().
	 */
        virtual int sync()
        {
	  //Wrap the text to the correct length
	  std::string ostring = string::linewrap<true>(str(), _linelength);

	  bool endl = *(ostring.end()-1) == '\n';
	  
	  //We strip the end "\n" for the search replace of the prefix
	  if (endl) ostring.erase(ostring.end() -1);
	  //Add the prefix to every newline
	  ostring = string::search_replace(ostring, "\n", "\n" + _prefix);
	  //Readd the final endline
	  if (endl) ostring += "\n";

	  //Write the result out
	  _output << _prefix + ostring;
	  _output.flush();
	  str("");
	  return 0;
        }
      };

      FormatingStreamBuf _buffer;

    public:
      /*! \brief Constructor.
       *
       * Warning: This class stores a pointer to the underlying ostream,
       * therefore the ostream must not fall out of scope before the
       * FormattedOStream does.
       *
       * \param ostream The underlying output stream or final destination of the formatted output.
       * \param prefix The string to replace all newline characters with.
       */
      inline FormattedOStream(const std::string & prefix = "",
			      std::ostream& ostream = std::cout,
			      const size_t linelength = 80):
	std::ostream(&_buffer),
	_buffer(prefix, ostream, linelength)
      {}
    
      /*! \brief Copy constructor.
       *
       * Typically std::ostream classes are not allowed copy
       * constructors, however as this is a ostream wrapper (and its
       * needed in DynamO) we allow it here.
       */
      FormattedOStream(const FormattedOStream& os):
	std::ostream(&_buffer),
	_buffer(os._buffer)
      {}
    };
  }
}
