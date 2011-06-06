/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
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

#include "1RSingle.hpp"
#include "../../simulation/particle.hpp"
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

CRSingle::CRSingle(const magnet::xml::Node& XML) 
{ operator<<(XML); }

bool 
CRSingle::isInRange(const Particle &part) const
{
  if (part.getID() == ID)
    return true;
  
    return false;
}

//The data output classes
void 
CRSingle::operator<<(const magnet::xml::Node& XML)
{
  if (strcmp(XML.getAttribute("Range"),"Single"))
    M_throw() << "Attempting to load CRSingle from non single";
    try {
      ID = XML.getAttribute("ID").as<size_t>();
    }
    catch (boost::bad_lexical_cast &)
      { M_throw() << "Failed a lexical cast in CRRange"; }
}

void 
CRSingle::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Range") << "Single"
      << xml::attr("ID") << ID;
}
