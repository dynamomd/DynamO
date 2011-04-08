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

#include "include.hpp"
#include "../../base/is_simdata.hpp"
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/exception.hpp>

xml::XmlStream& operator<<(xml::XmlStream& XML, 
			    const Units& g)
{
  g.outputXML(XML);
  return XML;
}

double 
Units::simVolume() const
{ 
  double vol = 1.0;
  for (size_t iDim = 0; iDim < NDIM; iDim++)
    vol *= Sim->aspectRatio[iDim];
  
  return vol;
}

Units* 
Units::getClass(const magnet::xml::Node& XML, const DYNAMO::SimData* Sim)
{
  if (!strcmp(XML.getAttribute("Type"),"Elastic")
      || !strcmp(XML.getAttribute("Type"),"HardSphere"))
    return new UHardSphere(XML, Sim);
  else if (!strcmp(XML.getAttribute("Type"),"Shear"))
    return new UShear(XML, Sim);
  else if (!strcmp(XML.getAttribute("Type"),"SW"))
    return new USquareWell(XML, Sim);
  else
   M_throw() << XML.getAttribute("Type")
	     << ", Unknown unit type";
}
