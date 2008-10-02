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

#ifndef __XMLWRITER_H_2CC1D410_9DE8_4452_B8F0_F987EF152E06
#define __XMLWRITER_H_2CC1D410_9DE8_4452_B8F0_F987EF152E06

#define XML_SPACING "  "

#define DYNAMO_XMLWrite

// disable terrible MSVC warnings which are issued when using STL
#ifdef	_MSC_VER
#pragma warning( disable : 4786 ) 
#pragma warning( disable : 4514 )
#endif

#include	<stack>
#include	<string>
#include	<sstream>

namespace xmlw {
  
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
	whatProlog, 
	whatTag, 
	whatTagEnd, 
	whatAttribute, 
	whatCharData
      }	what_type;
      
      what_type	what;
      std::string str;
      
      inline Controller(const Controller& c) : what(c.what), str(c.str) {}
      inline Controller(const what_type _what) : what(_what){}
      
      // use template constructor because string field <str> may be initialized 
      // from different sources: char*, std::string etc
      template<class t>
      inline Controller(const what_type _what, const t& _str): 
	what(_what), str(_str) {}
    };
    
    // XmlStream refers std::ostream object to perform actual output operations
    inline XmlStream(std::ostream& _s):
      state(stateNone), s(_s), prologWritten(false) {}
    
    inline XmlStream(const XmlStream &XML):
      state(XML.state), s(XML.s), prologWritten(XML.prologWritten) {}
    
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
    return XmlStream::Controller(XmlStream::Controller::whatProlog);
  }
  
  inline const XmlStream::Controller tag() {
    return XmlStream::Controller(XmlStream::Controller::whatTag);
  }
  
  inline const XmlStream::Controller tag(const char* const tag_name) {
    return XmlStream::Controller(XmlStream::Controller::whatTag, tag_name);
  }
  
  inline const XmlStream::Controller endtag() {
    return XmlStream::Controller(XmlStream::Controller::whatTagEnd);
  }
  
  inline const XmlStream::Controller endtag(const char* const tag_name) {
    return XmlStream::Controller(XmlStream::Controller::whatTagEnd, tag_name);
  }
  
  inline const XmlStream::Controller attr(const char* const attr_name) {
    return XmlStream::Controller(XmlStream::Controller::whatAttribute, attr_name);
  }
  
  inline const XmlStream::Controller chardata() {
    return XmlStream::Controller(XmlStream::Controller::whatCharData);
  }
  
}	// namespace


#endif  //  __XMLWRITER_H_2CC1D410_9DE8_4452_B8F0_F987EF152E06

