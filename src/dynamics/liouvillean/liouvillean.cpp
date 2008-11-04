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
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../../base/is_exception.hpp"
#include "../../base/is_simdata.hpp"
#include "../2particleEventData.hpp"
#include <boost/foreach.hpp>

xmlw::XmlStream& operator<<(xmlw::XmlStream& XML, const CLiouvillean& g)
{
  g.outputXML(XML);
  return XML;
}

CLiouvillean* 
CLiouvillean::loadClass(const XMLNode& XML, DYNAMO::SimData* tmp)
{
  if (!strcmp(XML.getAttribute("Type"),"Newtonian"))
    return new CLNewton(tmp);
  else 
    D_throw() << "Unknown type of Liouvillean encountered";
}

C2ParticleData 
CLiouvillean::runLineLineCollision() const
{ D_throw() << "Not implemented for this Liouvillean."; }

Iflt 
CLiouvillean::getLineLineCollision() const 
{ D_throw() << "Not implemented for this Liouvillean."; }
