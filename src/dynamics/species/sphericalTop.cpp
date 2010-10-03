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

#include "sphericalTop.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../units/units.hpp"
#include "../../base/is_simdata.hpp"

SpSphericalTop::SpSphericalTop(DYNAMO::SimData* tmp, CRange* nr, double nMass, 
			       std::string nName, unsigned int nID, double inertiaConst,
			       std::string nIName):
  SpInertia(tmp, "Species", IC_blue, nr, nMass, nName, nID, nIName),
  inertiaConstant(inertiaConst)
{}

SpSphericalTop::SpSphericalTop(const XMLNode& XML, DYNAMO::SimData* tmp, unsigned int nID):
  SpInertia(XML, tmp, nID)
{ operator<<(XML); }


void 
SpSphericalTop::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("InertiaConstant") 
      << inertiaConstant / Sim->dynamics.units().unitArea()
      << xml::attr("Mass") << mass / Sim->dynamics.units().unitMass()
      << xml::attr("Name") << spName
      << xml::attr("IntName") << intName
      << xml::attr("Type") << "SphericalTop"
      << range;
}

void 
SpSphericalTop::operator<<(const XMLNode& XML)
{
  Species::operator<<(XML);

  try {
    inertiaConstant 
      = boost::lexical_cast<double>(XML.getAttribute("InertiaConstant"))
      * Sim->dynamics.units().unitArea();
  } 
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in CSSphericalTop";
    }
}
