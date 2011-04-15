/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
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
#include <stack>
#include <string>
#include <sstream>

namespace xml {
  class XmlStream {
  public:
    // XML version constants
    static const int versionMajor = 1;
    static const int versionMinor = 0;
    
    // Internal helper class
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
    
    // XmlStream refers std::ostream object to perform actual output operations
    inline XmlStream(std::ostream& _s):
      state(stateNone), s(_s), prologWritten(false), FormatXML(false) {}
    
    inline XmlStream(const XmlStream &XML):
      state(XML.state), s(XML.s), prologWritten(XML.prologWritten), FormatXML(XML.FormatXML) {}
    
    inline ~XmlStream()
    {
      if (stateTagName == state) {
	s << "/>";
	state = stateNone;
      }
      while (tags.size())
	endTag(tags.top());
    }

    
    // this is the main working horse
    inline XmlStream& operator<<(const Controller& controller)
    {
      switch (controller._type) {
      case Controller::Prolog:
	if (!prologWritten && stateNone == state) {
	  s << "<?xml version=\"" << versionMajor << '.' << versionMinor << "\"?>\n";
	  prologWritten = true;
	}
	break;	//	Controller::whatProlog
	
      case Controller::Tag:
	closeTagStart();
	s << '<';
	if (controller.str.empty()) {
	  clearTagName();
	  state = stateTagName;
	}
	else {
	  s << controller.str;
	  tags.push(controller.str);
	  state = stateTag;
	}
	break;	//	Controller::whatTag
      
      case Controller::TagEnd:
	endTag(controller.str);
	break;	//	Controller::whatTagEnd
      
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

    // default behaviour - delegate object output to std::stream
    template<class t>
    XmlStream& operator<<(const t& value) {
      if (stateTagName == state)
	tagName << value;
      s << value;
      return *this;
    }

    inline std::ostream& getUnderlyingStream() { return s; }

    inline void setFormatXML(const bool& tf) { FormatXML = tf; }
    
  private:
    // state of the stream 
    typedef enum 
      {
	stateNone, 
	stateTag, 
	stateAttribute, 
	stateTagName, 
	stateCharData
      }	state_type;
    
    // tag name stack
    typedef std::stack<std::string>	tag_stack_type;
    
    tag_stack_type	tags;
    state_type	state;
    std::ostream&	s;
    bool	prologWritten;
    std::ostringstream	tagName;
    bool        FormatXML;
    
    inline void clearTagName() {
      tagName.str(std::string());
    }
    
    // Close current tag
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

    // Close tag (may be with closing all of its children)
    inline void endTag(const std::string& tag)
    {
      bool brk = false;
    
      while (tags.size() > 0 && !brk) {
	if (stateNone == state)
	  s << "</" << tags.top() << ">\n";
	else {
	  closeTagStart(true);
	  state = stateNone;
	}
	brk = tag.empty() || tag == tags.top();
	tags.pop();
      }
    }    
  };	//	class XmlStream
  
  // Helper functions, they may be simply overwritten
  // E.g. you may use std::string instead of const char*
  
  inline const XmlStream::Controller prolog() {
    return XmlStream::Controller(XmlStream::Controller::Prolog);
  }
  
  inline const XmlStream::Controller tag() {
    return XmlStream::Controller(XmlStream::Controller::Tag);
  }
  
  inline const XmlStream::Controller tag(const std::string& tag_name) {
    return XmlStream::Controller(XmlStream::Controller::Tag, tag_name);
  }
  
  inline const XmlStream::Controller endtag() {
    return XmlStream::Controller(XmlStream::Controller::TagEnd);
  }
  
  inline const XmlStream::Controller endtag(const std::string& tag_name) {
    return XmlStream::Controller(XmlStream::Controller::TagEnd, tag_name);
  }
  
  inline const XmlStream::Controller attr(const std::string& attr_name) {
    return XmlStream::Controller(XmlStream::Controller::Attribute, attr_name);
  }
  
  inline const XmlStream::Controller chardata() {
    return XmlStream::Controller(XmlStream::Controller::CharData);
  }
}


