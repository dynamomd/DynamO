/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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

#include <dynamo/dynamics/interactions/nullInteraction.hpp>
#include <dynamo/dynamics/interactions/intEvent.hpp>
#include <dynamo/dynamics/2particleEventData.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <cstring>

namespace dynamo {
  INull::INull(dynamo::SimData* tmp, C2Range* nR, std::string name):
    Interaction(tmp, nR) { intName = name; }

  INull::INull(const magnet::xml::Node& XML, dynamo::SimData* tmp):
    Interaction(tmp, NULL)
  {
    operator<<(XML);
  }

  void 
  INull::initialise(size_t nID)
  { ID=nID; }

  void 
  INull::operator<<(const magnet::xml::Node& XML)
  { 
    if (std::strcmp(XML.getAttribute("Type"),"Null"))
      M_throw() << "Attempting to load NullInteraction from " 
		<< XML.getAttribute("Type") <<" entry";
  
    Interaction::operator<<(XML);
  
    try 
      { intName = XML.getAttribute("Name"); }
    catch (boost::bad_lexical_cast &)
      {
	M_throw() << "Failed a lexical cast in CINull";
      }
  }

  IntEvent 
  INull::getEvent(const Particle &p1, const Particle &p2) const 
  { 
    return IntEvent(p1, p2, HUGE_VAL, NONE, *this);
  }

  void
  INull::runEvent(Particle&, Particle&, const IntEvent&) const
  { 
    M_throw() << "Null event trying to run a collision!"; 
  }
   
  void 
  INull::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "Null"
	<< magnet::xml::attr("Name") << intName
	<< *range;
  }
}

