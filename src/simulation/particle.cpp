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

#include "particle.hpp"
#include "../base/is_exception.hpp"
#include "../extcode/xmlwriter.hpp"
#include "../extcode/xmlParser.h"
#include "../datatypes/vector.xml.hpp"
#include <ostream>

CParticle::CParticle(const XMLNode& XML, unsigned long nID):
  ID(nID),
  pecTime(0.0)
{
  XMLNode xBrowseNode = XML.getChildNode("P");
  posVector << xBrowseNode;
  
  xBrowseNode = XML.getChildNode("V");
  velVector << xBrowseNode;
}

xmlw::XmlStream& operator<<(xmlw::XmlStream& XML, 
			    const CParticle& particle)
{
  XML << xmlw::attr("ID") << particle.ID
      << xmlw::tag("P")
      << (particle.posVector)
      << xmlw::endtag("P")
      << xmlw::tag("V")
      << (particle.velVector)
      << xmlw::endtag("V");
  
  return XML;
}
