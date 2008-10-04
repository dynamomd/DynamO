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

#include "include.hpp"
#include "../../extcode/xmlParser.h"
#include "../../extcode/xmlwriter.hpp"
#include "../../base/is_exception.hpp"

xmlw::XmlStream& operator<<(xmlw::XmlStream& XML, 
			    const CBC& g)
{
  g.outputXML(XML);
  return XML;
}


CBC* 
CBC::loadClass(const XMLNode &XML, DYNAMO::SimData* tmp)
{
  if (!strcmp(XML.getAttribute("Boundary"),"Null"))
    return new CNullBC(tmp);
  else if (!strcmp(XML.getAttribute("Shape"),"Square"))
    {
      if (!strcmp(XML.getAttribute("Boundary"),"PBC"))
	return new CSPBC(tmp);
      else if (!strcmp(XML.getAttribute("Boundary"),"LE"))
	return new CSLEBC(XML,tmp);
      else 
	D_throw() << "Unknown type of square boundary encountered";
    } 
  else if (!strcmp(XML.getAttribute("Shape"),"Rectangular"))
    {
      if (!strcmp(XML.getAttribute("Boundary"),"PBC"))
	return new CRPBC(tmp);
      else if (!strcmp(XML.getAttribute("Boundary"),"NoXPBC"))
	return new CRNoXPBC(tmp);
      else if (!strcmp(XML.getAttribute("Boundary"),"LE"))
	return new CRLEBC(XML,tmp);
      else 
	D_throw() << "Unknown type of rectangular boundary encountered";
    }
  else
    D_throw() << "Unknown shape of Boundary encountered";

}
