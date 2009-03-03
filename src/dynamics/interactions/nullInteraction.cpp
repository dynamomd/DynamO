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

#include "nullInteraction.hpp"
#include "../../base/is_exception.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../../dynamics/interactions/intEvent.hpp"
#include "../2particleEventData.hpp"
#include <cstring>

CINull::CINull(DYNAMO::SimData* tmp, C2Range* nR):
  CInteraction(tmp, nR) {}

CINull::CINull(const XMLNode& XML, DYNAMO::SimData* tmp):
  CInteraction(tmp,NULL)
{
  operator<<(XML);
}

void 
CINull::initialise(size_t nID)
{ ID=nID; }

void 
CINull::operator<<(const XMLNode& XML)
{ 
  if (std::strcmp(XML.getAttribute("Type"),"Null"))
    D_throw() << "Attempting to load NullInteraction from " 
	      << XML.getAttribute("Type") <<" entry";
  
  range.set_ptr(C2Range::loadClass(XML,Sim));
  
  try 
    {
      intName = XML.getAttribute("Name");
    }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in CINull";
    }
}

Iflt 
CINull::maxIntDist() const 
{ return 0; }

Iflt 
CINull::hardCoreDiam() const 
{ return 0; }

void 
CINull::rescaleLengths(Iflt) 
{}

CInteraction* 
CINull::Clone() const 
{ return new CINull(*this); }
  
CIntEvent 
CINull::getEvent(const CParticle &p1, const CParticle &p2) const 
{ 
  return CIntEvent(p1, p2, HUGE_VAL, NONE, *this);
}

void
CINull::runEvent(const CParticle&, const CParticle&, const CIntEvent&) const
{ 
  D_throw() << "Null event trying to run a collision!"; 
}
   
void 
CINull::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "Null"
      << xmlw::attr("Name") << intName
      << range;
}

void
CINull::checkOverlaps(const CParticle&, const CParticle&) const
{}
