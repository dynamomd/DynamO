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

#include "xmlwriter.hpp"

namespace xmlw
{
  // Before destroying check whether all the open tags are closed
  XmlStream::~XmlStream() 
  {
    if (stateTagName == state) {
      s << "/>";
      state = stateNone;
    }
    while (tags.size())
      endTag(tags.top());
  }
  

  XmlStream& 
  XmlStream::operator<<(const Controller& controller) 
  {
    switch (controller.what) {
    case Controller::whatProlog:
      if (!prologWritten && stateNone == state) {
	  s << "<?xml version=\"" << versionMajor << '.' << versionMinor << "\"?>\n";
	  prologWritten = true;
      }
      break;	//	Controller::whatProlog
      
    case Controller::whatTag:
      closeTagStart();
      if (FormatXML) for (unsigned int i = 0; i < tags.size(); i++) s << XML_SPACING;
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
      
    case Controller::whatTagEnd:
      endTag(controller.str);
      break;	//	Controller::whatTagEnd
      
    case Controller::whatAttribute:
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
      // else throw some error - unexpected attribute (out of any tag)
      
      break;	//	Controller::whatAttribute
      
    case Controller::whatCharData:
      closeTagStart();
      state = stateNone;
      break;	//	Controller::whatCharData
    }
    
    return	*this;
  }

  // Close current tag
  void 
  XmlStream::closeTagStart(bool self_closed) 
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

  void 
  XmlStream::endTag(const std::string& tag) 
  {
    bool brk = false;
    
    while (tags.size() > 0 && !brk) {
      if (stateNone == state)
	{
	  if (FormatXML) 
	    for (unsigned int i = 0; i < tags.size()-1; i++) 
	      s << XML_SPACING;

	  s << "</" << tags.top() << '>' << '\n';
	}
      else {
	closeTagStart(true);
	state = stateNone;
      }
      brk = tag.empty() || tag == tags.top();
      tags.pop();
    }
  }
  
}
