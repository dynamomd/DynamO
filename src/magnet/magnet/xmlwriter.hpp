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
#include <tr1/memory>
#include <stack>
#include <string>
#include <sstream>

namespace magnet {
  namespace xml {
    /*! \brief A class which behaves like an output stream for XML output.
     */
    class XmlStream {
    public:
      //! \brief Major XML version constant.
      static const int versionMajor = 1;
      //! \brief Minor XML version constant.
      static const int versionMinor = 0;
    
      //! \brief Internal type used to modify the state of the XML stream.
      struct Controller {
	typedef enum {
	  Prolog, 
	  Tag, 
	  TagEnd, 
	  Attribute, 
	  CharData
	}	ControllerType;
      
	ControllerType _type;
	std::string str;
      
	inline Controller(const Controller& c) : _type(c._type), str(c.str) {}
	inline Controller(const ControllerType type) : _type(type){}
	inline Controller(const ControllerType type, const std::string& _str):
	  _type(type), str(_str) {}
      };
    
      //! \brief Constructs an XmlStream from a std::ostream object.
      inline XmlStream(std::ostream& _s):
	state(stateNone), s(_s), prologWritten(false), FormatXML(false) {}
    
      //! \brief Copy constructor.
      inline XmlStream(const XmlStream &XML):
	state(XML.state), s(XML.s), prologWritten(XML.prologWritten), FormatXML(XML.FormatXML) {}
    
      //! \brief Destructor.
      inline ~XmlStream()
      {
	if (stateTagName == state) {
	  s << "/>";
	  state = stateNone;
	}
	while (tags.size())
	  endTag(tags.top());
      }

    
      /*! \brief Main insertion operator which changes the state of the
       * XmlStream.
       */
      inline XmlStream& operator<<(const Controller& controller)
      {
	switch (controller._type) {
	case Controller::Prolog:
	  if (!prologWritten && stateNone == state) {
	    s << "<?xml version=\"" << versionMajor << '.' << versionMinor << "\"?>\n";
	    prologWritten = true;
	  }
	  break;
	case Controller::Tag:
	  closeTagStart();
	  for (size_t i(0); i < tags.size(); ++i)
	    s << "  ";
	  s << '<' << controller.str;
	  tags.push(controller.str);
	  state = stateTag;
	  break;
	case Controller::TagEnd:
	  endTag(controller.str);
	  break;
	case Controller::Attribute:
	  switch (state) {
	  case stateTagName:
	    tags.push(tagName.str());
	    break;
	  case stateAttribute:
	    s << '\"';
	  default:
	    break;
	  }
      
	  if (stateNone != state) {
	    s << ' ' << controller.str << "=\"";
	    state = stateAttribute;
	  }
	  break;//Controller::whatAttribute
	case Controller::CharData:
	  closeTagStart();
	  state = stateNone;
	  break;//Controller::whatCharData
	}
	return	*this;
      }

      /*! \brief Default insertion operator, just delegates the passed
       * object to the underlying std::stream.
       */
      template<class T>
      XmlStream& operator<<(const T& value) {
	if (stateTagName == state)
	  tagName << value;
	s << value;
	return *this;
      }

      /*! \brief Specialisation for pointers.
	
	Calls the output operator on the dereferenced object
       */
      template<class T>
      XmlStream& operator<<(const std::tr1::shared_ptr<T>& value) {
	if (stateTagName == state)
	  tagName << value;
	return ((*this) << (*value));
      }

      //! \brief Returns the underlying output stream.
      inline std::ostream& getUnderlyingStream() { return s; }

      /*! \brief Enables or disables automatic formatting of the
       * outputted XML.
       */
      inline void setFormatXML(const bool& tf) { FormatXML = tf; }
    
    private:
      //! \brief Enum types used to track the current state of the XmlStream.
      typedef enum 
	{
	  stateNone, 
	  stateTag, 
	  stateAttribute, 
	  stateTagName, 
	  stateCharData
	}	state_type;
    
      //! \brief Stack of parent XML nodes above the current node.
      typedef std::stack<std::string>	tag_stack_type;
    
      tag_stack_type	tags;
      state_type	state;
      std::ostream&	s;
      bool	prologWritten;
      std::ostringstream	tagName;
      bool        FormatXML;
    
      //! \brief Closes the current tag.
      inline void closeTagStart(bool self_closed = false)
      {
	if (stateTagName == state)
	  tags.push(tagName.str());
      
	// note: absence of 'break's is not an error
	switch (state) {
	case stateAttribute:
	  s << '\"';
	case stateTagName:
	case stateTag:
	  if (self_closed)
	    s << '/';
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
	  if (stateNone == state)
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
      }    
    };
  
    /*! \brief Stream manipulator causing the XmlStream to output the
     * XML prolog.
     */
    inline const XmlStream::Controller prolog() { return
	XmlStream::Controller(XmlStream::Controller::Prolog); }
  
    /*! \brief Stream manipulator to add an XML tag to the XmlStream.
     */
    inline const XmlStream::Controller tag(const std::string& tag_name) {
      return XmlStream::Controller(XmlStream::Controller::Tag, tag_name);
    }
  
    /*! \brief Stream manipulator to close an XML tag in the XmlStream.
     *
     * This will close all XML tags until it closes one tag with the
     * passed name.
     */
    inline const XmlStream::Controller endtag(const std::string& tag_name) {
      return XmlStream::Controller(XmlStream::Controller::TagEnd, tag_name);
    }
  
    /*! \brief Stream manipulator to add an attribute with the passed
     * name to the XmlStream.
     *
     * This manipulates the stream state such that the next passed value
     * is used as the value of the attribute.
     */
    inline const XmlStream::Controller attr(const std::string& attr_name) {
      return XmlStream::Controller(XmlStream::Controller::Attribute, attr_name);
    }
  
    /*! \brief Stream manipulator to switch the stream state to output
     * all further inserted values to the contents of the current XML
     * tag.
     *
     * This mode is ended by the next \ref endtag stream modifier.
     */
    inline const XmlStream::Controller chardata() {
      return XmlStream::Controller(XmlStream::Controller::CharData);
    }
  }
}
