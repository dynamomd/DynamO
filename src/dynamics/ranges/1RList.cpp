/*  DYNAMO:- Event driven molecular dynamics simulator 
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
#include "../../extcode/xmlParser.h"
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <magnet/exception.hpp>
#include <magnet/xmlwriter.hpp>

CRList::CRList(const XMLNode& XML) 
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
CRList::operator<<(const XMLNode& XML)
{
  if (strcmp(XML.getAttribute("Range"),"List"))
    M_throw() << "Attempting to load CRList from non list";
  try {
    
    XMLNode xSubNode;
    for (long i=0; i < XML.nChildNode("ID"); i++)
      IDs.push_back(boost::lexical_cast<unsigned long>
		    (XML.getChildNode("ID",i).getAttribute("val")));
  }
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in CRList";
    }
}

void 
CRList::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Range") << "List";
  BOOST_FOREACH(unsigned long ID, IDs)
    XML << xml::tag("ID") << xml::attr("val") << ID << xml::endtag("ID");
}
