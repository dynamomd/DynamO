/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>
    Copyright (C) 2003  Oblutus <http://iridia.ulb.ac.be/~fvandenb/tools/xmlParser.html>

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
#include <memory>
#include <stack>
#include <string>
#include <sstream>
#include <magnet/exception.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/chain.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/copy.hpp>

namespace magnet {
  namespace xml {
    namespace detail {
      /*! \brief  A custom file sink for the boost iostreams library.
	
	Boost Iostreams library does not correctly check if files are
	in a good state for writing, during the write itself. Thus, if
	the disk becomes full, or the program does not have write
	permission, the writing will hang in an infinite loop! This
	file sink fixes that.
      */
      struct safe_ofstream_sink {
	typedef char char_type;
	typedef boost::iostreams::sink_tag category;
	
	safe_ofstream_sink(std::string filename):
	  _ofs(new std::ofstream(filename)),
	  _filename(filename)
	{
	  if (not (*_ofs))
	    M_throw() << "Failed to open " << _filename << " for writing.";
	}

	~safe_ofstream_sink() {
	  if (_ofs and _ofs->is_open()) {
	    _ofs->flush();
	    if (not (*_ofs))
	      M_throw() << "Failed to flush " << _filename << ".";
	  }
	}

	std::streamsize write(const char* s, std::streamsize n)
	{
	  if (not _ofs)
	    M_throw() << "File handle lost during copying!";
	  
	  _ofs->write(s, n);
	  if (not (*_ofs))
            M_throw() << "Failed while writing to" << _filename;

	  return n;
	}
	
	std::shared_ptr<std::ofstream> _ofs;
	std::string _filename;
      };

    }
		    
    
    namespace io = boost::iostreams;

    /*! \brief A class which behaves like an output stream for XML output.
     */
    class XmlStream {
    public:

      //! \brief Internal type used to modify the state of the XML stream.
      struct Controller {
	typedef enum {
	  Prolog, 
	  Tag, 
	  TagEnd, 
	  Attribute, 
	  CharData
	} ControllerType;
	ControllerType _type;
	std::string _str;
      
	inline Controller(const ControllerType type) : _type(type){}
	inline Controller(const ControllerType type, const std::string& str):
	  _type(type), _str(str) {}
      };
    
      inline XmlStream(std::string filename):
	state(stateNone), prologWritten(false), FormatXML(false)
      {
	if (std::string(filename.end() - 4, filename.end()) == ".bz2") {
#ifdef DYNAMO_bzip2_support
	  s.push(io::bzip2_compressor());
#else
	  M_throw() << "bz2 compressed file support was not built in! (only available on linux)";
#endif
	}
	s.push(detail::safe_ofstream_sink(filename));
      }
        
      inline ~XmlStream() { while (tags.size()) endTag(tags.top()); }

      /*! \brief Main insertion operator which changes the state of
        the XmlStream.
       */
      inline XmlStream& operator<<(const Controller& controller)
      {
	switch (controller._type) {
	case Controller::Prolog:
	  if (prologWritten) M_throw() << "XML prolog already written.";
	  if (state != stateNone) M_throw() << "Incorrect state to write XML prolog.";
	  s << "<?xml version=\"1.0\"?>\n";
	  prologWritten = true;
	  break;
	case Controller::Tag:
	  closeTagStart();
	  for (size_t i(0); i < tags.size(); ++i) s << "  ";
	  s << '<' << controller._str;
	  tags.push(controller._str);
	  state = stateTag;
	  break;
	case Controller::TagEnd:
	  endTag(controller._str);
	  break;
	case Controller::Attribute:
	  if (state == stateAttribute) s << '\"'; //Close open attribute
	  if (state != stateNone) {
	    s << ' ' << controller._str << "=\"";
	    state = stateAttribute;
	  }
	  break;
	case Controller::CharData:
	  closeTagStart();
	  state = stateNone;
	  break;
	}
	return	*this;
      }

      /*! \brief Default insertion operator, just delegates the passed
	object to the underlying std::stream.
       */
      template<class T>
      friend XmlStream& operator<<(XmlStream& XML, const T& value) {
	static_cast<std::ostream &>(XML.s) << value;
	return XML;
      }

      /*! \brief Specialisation for pointers. */
      template<class T>
      friend XmlStream& operator<<(XmlStream& XML, const std::shared_ptr<T>& value) {
	return XML << (*value);
      }

      //! \brief Returns the underlying output stream.
      inline std::ostream& getUnderlyingStream() { return s; }

      /*! \brief Enables or disables automatic formatting of the
        outputted XML.
       */
      inline void setFormatXML(const bool& tf) { FormatXML = tf; }
    
    private:
      //! \brief Enum types used to track the current state of the XmlStream.
      typedef enum 
	{
	  stateNone, 
	  stateTag, 
	  stateAttribute, 
	  stateCharData
	}	state_type;
    
      //! \brief Stack of parent XML nodes above the current node.
      typedef std::stack<std::string>	tag_stack_type;
    
      tag_stack_type	tags;
      state_type	state;
      io::filtering_ostream s;
      bool	prologWritten;
      std::ostringstream	tagName;
      bool        FormatXML;
    
      //! \brief Closes the current tag.
      inline void closeTagStart(bool self_closed = false)
      {
	// note: absence of 'break's is not an error
	switch (state) {
	case stateAttribute:
	  s << '\"';
	case stateTag:
	  if (self_closed) s << '/';
	  s << '>' << '\n';
	default:
	  break;
	}
      }

      //! \brief Closes the named tag and may close all of its children.
      inline void endTag(const std::string& tag)
      {
	bool brk = false;
    
	while (tags.size() > 0 && !brk) {
	  if (state == stateNone)
	    {
	      for (size_t i(1); i < tags.size(); ++i)
		s << "  ";
	      s << "</" << tags.top() << ">\n";
	    }
	  else {
	    closeTagStart(true);
	    state = stateNone;
	  }
	  brk = tag.empty() || tag == tags.top();
	  tags.pop();
	}

	//Check if we ran out of tags
	if (tags.empty() && !brk)
	  M_throw() << "Could not find tag \"" << tag << "\" to close!";
      }
    };
  
    /*! \brief Stream manipulator causing the XmlStream to output the
      XML prolog.
     */
    inline const XmlStream::Controller prolog() { return
	XmlStream::Controller(XmlStream::Controller::Prolog); }
  
    /*! \brief Stream manipulator to add an XML tag to the XmlStream.
     */
    inline const XmlStream::Controller tag(const std::string& tag_name) {
      return XmlStream::Controller(XmlStream::Controller::Tag, tag_name);
    }
  
    /*! \brief Stream manipulator to close an XML tag in the XmlStream.
     
      This will close all XML tags until it closes one tag with the
      passed name.
     */
    inline const XmlStream::Controller endtag(const std::string& tag_name) {
      return XmlStream::Controller(XmlStream::Controller::TagEnd, tag_name);
    }
  
    /*! \brief Stream manipulator to add an attribute with the passed
      name to the XmlStream.
     
      This manipulates the stream state such that the next passed value
      is used as the value of the attribute.
     */
    inline const XmlStream::Controller attr(const std::string& attr_name) {
      return XmlStream::Controller(XmlStream::Controller::Attribute, attr_name);
    }
  
    /*! \brief Stream manipulator to switch the stream state to output
      all further inserted values to the contents of the current XML
      tag.
     
      This mode is ended by the next \ref endtag stream modifier.
     */
    inline const XmlStream::Controller chardata() {
      return XmlStream::Controller(XmlStream::Controller::CharData);
    }
  }
}
