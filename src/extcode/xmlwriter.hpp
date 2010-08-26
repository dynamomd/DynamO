/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

//////////////////////////////////////////////////////////////////////
//
// xmlparser.h: interface & implementation for the XmlStream class.
// 
// Author: Oboltus, December 2003
// Formatted output by Marcus Bannerman (http://marcusbannerman.co.uk)
//
// http://iridia.ulb.ac.be/~fvandenb/tools/xmlParser.html
//
// This code is provided "as is", with absolutely no warranty expressed
// or implied. Any use is at your own risk.
//
//////////////////////////////////////////////////////////////////////

#pragma once

#include	<stack>
#include	<string>
#include	<sstream>

namespace xml {
  //! The spacing unit to use when a tag is opened
  const std::string XML_SPACING("  ");
  
  class XmlStream {
  public:
    // XML version constants
    enum {
      versionMajor = 1, 
      versionMinor = 0
    };
    
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
      
      // use template constructor because string field <str> may be initialized 
      // from different sources: char*, std::string etc
      template<class t>
      inline Controller(const ControllerType type, const t& _str):
	_type(type), str(_str) {}
    };
    
    // XmlStream refers std::ostream object to perform actual output operations
    inline XmlStream(std::ostream& _s):
      state(stateNone), s(_s), prologWritten(false), FormatXML(false) {}
    
    inline XmlStream(const XmlStream &XML):
      state(XML.state), s(XML.s), prologWritten(XML.prologWritten), FormatXML(XML.FormatXML) {}
    
    ~XmlStream();
    
    // this is the main working horse
    XmlStream& operator<<(const Controller&);
    
    // default behaviour - delegate object output to std::stream
    template<class t>
    XmlStream& operator<<(const t& value) {
      if (stateTagName == state)
	tagName << value;
      s << value;
      return *this;
    }

    std::ostream& getUnderlyingStream() { return s; }

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
    
    // I don't know any way easier (legal) to clear std::stringstream...
    inline void clearTagName() {
      const std::string	empty_str;
      tagName.rdbuf()->str(empty_str);
    }
    
    // Close current tag
    void closeTagStart(bool self_closed = false);
    
    // Close tag (may be with closing all of its children)
    void endTag(const std::string&);

  };	//	class XmlStream
  
  // Helper functions, they may be simply overwritten
  // E.g. you may use std::string instead of const char*
  
  inline const XmlStream::Controller prolog() {
    return XmlStream::Controller(XmlStream::Controller::Prolog);
  }
  
  inline const XmlStream::Controller tag() {
    return XmlStream::Controller(XmlStream::Controller::Tag);
  }
  
  inline const XmlStream::Controller tag(const char* const tag_name) {
    return XmlStream::Controller(XmlStream::Controller::Tag, tag_name);
  }
  
  inline const XmlStream::Controller endtag() {
    return XmlStream::Controller(XmlStream::Controller::TagEnd);
  }
  
  inline const XmlStream::Controller endtag(const char* const tag_name) {
    return XmlStream::Controller(XmlStream::Controller::TagEnd, tag_name);
  }
  
  inline const XmlStream::Controller attr(const char* const attr_name) {
    return XmlStream::Controller(XmlStream::Controller::Attribute, attr_name);
  }
  
  inline const XmlStream::Controller chardata() {
    return XmlStream::Controller(XmlStream::Controller::CharData);
  }
  
}	// namespace


