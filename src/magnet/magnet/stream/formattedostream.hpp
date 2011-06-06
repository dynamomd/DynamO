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
#include <magnet/stream/console_specials.hpp>
#include <iostream>
#include <sstream>

namespace magnet
{
  namespace stream {
    /*! \brief A class which is used to mark the end of a block of
     * formatted text.
     *
     * This class is used like a std::endl class, as it causes all
     * text in the buffer to be formatted and flushed to the
     * underlying _ostream.
     */
    struct EndBlock {};

    /*! \brief This class wraps an std::ostream to provide automatic
     * formatting.
     *
     * The primary purpose of this class is to provide formatted
     * output for classes. The idea is that the output of each class
     * might be prefixed with some identifying information, and long
     * lines will be automatically wrapped. To acheive this, the end
     * of every chunk of information passed to this class must finish
     * with a FormattedOStream::EndBlock. e.g.
     *
     * \code 
     * stream::FormattedOStream os(...);
     * os << "Some long text as part of a block of output, plus a number " << 20 
     *    << "but always finished with a " << stream::EndBlock; \endcode
     */
    class FormattedOStream
    {
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
      inline FormattedOStream(std::ostream& ostream = std::cout, 
			      const std::string & prefix = "", 
			      const size_t linelength = 80):
	_ostream(&ostream),
	_prefix(prefix),
	_linelength(linelength - prefix.size())
      {
	//Add a reset to make sure any prefix formatting is not carried
	//over into the string
	_prefix = "\n" + _prefix;
      }
    
      /*! \brief The main workhorse for the FormattedOStream operator.
       *
       * This stores all passed types in the _ostring
       * std::stringstream object. These strings will be formatted
       * once a std::endl object is recieved, by the std::endl
       * specialization.
       */
      template<class T>
      inline FormattedOStream& operator<<(T m)
      {
	_ostring << m;
	return *this;
      }

      /*! \brief Associates the Stream_Operator with a stream. 
       *
       * This function lets the formatting class be used like a stream
       * modifier.
       */
      inline friend FormattedOStream& operator<<(std::ostream &os, FormattedOStream& SO)
      {
	SO._ostream = &os;
	return SO;
      }

      /*! \brief Get the underlying std::ostream object. */
      inline std::ostream& getStream() const { return *_ostream; }
    
    protected:
      /*! \brief Name to insert after newlines.*/
      std::string _prefix;
    
      /*! \brief The underlying std::ostream to send the formatted
        strings to. */
      std::ostream *_ostream;

      /*! \brief The buffer, where all the unformatted strings are
       *  collected before a std::endl flushes them.
       */
      std::ostringstream _ostring;

      /*! \brief The maximum length of a formatted line before it is
       *   wrapped.
       */      
      size_t _linelength;
    };
    
    /*! \brief The actual output function.
     *
     * When a EndBlock is sent to the stream, we actually format and
     * output the text.
     */
    template<>
    inline FormattedOStream&
    FormattedOStream::operator<< <EndBlock> (EndBlock)
    {
      //Add the initial newline
      std::string output(std::string("\n") + _ostring.str());
      //Wrap the text to the correct length
      output = string::linewrap<true>(output, _linelength);
      //Add the prefix to every newline
      output = string::searchReplace(output, "\n", _prefix);
      //Write the result out
      *_ostream << output;
      return *this;
    }
  }
}
