/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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
#include "../../datatypes/vector.hpp"
#include <boost/lexical_cast.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../../base/is_exception.hpp"
#include "../../base/is_simdata.hpp"

xmlw::XmlStream& operator<<(xmlw::XmlStream& XML, 
			    const Units& g)
{
  g.outputXML(XML);
  return XML;
}

Iflt 
Units::simVolume() const
{ 
  Iflt vol = 1.0;
  for (size_t iDim = 0; iDim < NDIM; iDim++)
    vol *= Sim->aspectRatio[iDim];
  
  return vol;
}

Units* 
Units::loadUnits(const XMLNode &XML, const DYNAMO::SimData* Sim)
{
  if (!strcmp(XML.getAttribute("Type"),"Elastic")
      || !strcmp(XML.getAttribute("Type"),"HardSphere"))
    return new UHardSphere(XML, Sim);
  else if (!strcmp(XML.getAttribute("Type"),"Shear"))
    return new UShear(XML, Sim);
  else if (!strcmp(XML.getAttribute("Type"),"SW"))
    return new USquareWell(XML, Sim);
  else
   D_throw() << "Could not recognise the type of the units";
}
