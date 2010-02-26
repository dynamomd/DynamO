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

#include "1RRange.hpp"
#include "../../simulation/particle.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../../base/is_exception.hpp"

CRRange::CRRange(const XMLNode& XML) 
{ operator<<(XML); }

void 
CRRange::operator<<(const XMLNode& XML)
{
  if (strcmp(XML.getAttribute("Range"),"Ranged"))
    D_throw() << "Attempting to load CRRange from non range";
  try {
    startID = boost::lexical_cast<unsigned long>(XML.getAttribute("Start"));
    endID = boost::lexical_cast<unsigned long>(XML.getAttribute("End"));
  }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in CRRange";
    }
}

void 
CRRange::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Range") << "Ranged"
      << xmlw::attr("Start") << startID
      << xmlw::attr("End") << endID;
}
