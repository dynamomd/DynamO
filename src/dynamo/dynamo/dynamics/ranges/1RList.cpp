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

#include "1RList.hpp"
#include "../../simulation/particle.hpp"
#include <magnet/exception.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>

CRList::CRList(const magnet::xml::Node& XML) 
{ operator<<(XML); }

bool 
CRList::isInRange(const Particle &part) const
{
  BOOST_FOREACH(const unsigned long ID, IDs)
    if (part.getID() == ID)
      return true;
  return false;
}

void 
CRList::operator<<(const magnet::xml::Node& XML)
{
  if (strcmp(XML.getAttribute("Range"),"List"))
    M_throw() << "Attempting to load CRList from non list";
  try {
    
    for (magnet::xml::Node node = XML.getNode("ID"); node.valid(); ++node)
      IDs.push_back(node.getAttribute("val").as<unsigned long>());
  }
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in CRList";
    }
}

void 
CRList::outputXML(magnet::xml::XmlStream& XML) const
{
  XML << magnet::xml::attr("Range") << "List";
  BOOST_FOREACH(unsigned long ID, IDs)
    XML << magnet::xml::tag("ID") << magnet::xml::attr("val") << ID << magnet::xml::endtag("ID");
}
