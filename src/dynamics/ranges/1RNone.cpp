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

#include "1RNone.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"

CRNone::CRNone(const XMLNode& XML)
{ operator<<(XML); }

void 
CRNone::operator<<(const XMLNode& XML)
{
  if (strcmp(XML.getAttribute("Range"),"None"))
    M_throw() << "Attempting to load CRNone from non None type";
}

void 
CRNone::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Range") << "All";
}


