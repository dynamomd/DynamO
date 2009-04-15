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

#include "sphericalTop.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../units/units.hpp"
#include "../../base/is_simdata.hpp"

CSSphericalTop::CSSphericalTop(DYNAMO::SimData* tmp, CRange* nr, Iflt nMass, 
			       std::string nName, unsigned int nID, Iflt inertiaConst,
			       std::string nIName):
  CSpecInertia(tmp, "Species", IC_blue, nr, nMass, nName, nID, nIName),
  inertiaConstant(inertiaConst)
{}

CSSphericalTop::CSSphericalTop(const XMLNode& XML, DYNAMO::SimData* tmp, unsigned int nID):
  CSpecInertia(XML, tmp, nID)
{ operator<<(XML); }


void 
CSSphericalTop::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("InertiaConstant") 
      << inertiaConstant / Sim->Dynamics.units().unitArea()
      << xmlw::attr("Mass") << mass / Sim->Dynamics.units().unitMass()
      << xmlw::attr("Name") << spName
      << xmlw::attr("IntName") << intName
      << xmlw::attr("Type") << "SphericalTop"
      << range;
}

void 
CSSphericalTop::operator<<(const XMLNode& XML)
{
  CSpecies::operator<<(XML);

  try {
    inertiaConstant 
      = boost::lexical_cast<Iflt>(XML.getAttribute("InertiaConstant"))
      * Sim->Dynamics.units().unitArea();
  } 
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in CSSphericalTop";
    }
}
