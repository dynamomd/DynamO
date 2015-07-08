/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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

#include <dynamo/species/sphericalTop.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/BC/LEBC.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
  SpSphericalTop::SpSphericalTop(dynamo::Simulation* tmp, IDRange* nr, double nMass, 
				 std::string nName, unsigned int nID, double inertiaConst):
    SpInertia(tmp, nr, nMass, nName, nID),
    inertiaConstant(inertiaConst)
  {}

  SpSphericalTop::SpSphericalTop(const magnet::xml::Node& XML, dynamo::Simulation* Sim, 
				 unsigned int nID):
    SpInertia(XML, Sim, nID)
  { operator<<(XML); }


  void 
  SpSphericalTop::outputXML(magnet::xml::XmlStream& XML, std::string type) const
  {
    XML << magnet::xml::attr("InertiaConstant") 
	<< inertiaConstant / Sim->units.unitInertia()
	<< magnet::xml::attr("Mass") << _mass->getName()
	<< magnet::xml::attr("Name") << spName
	<< magnet::xml::attr("Type") << type
	<< range;
  }
  
  double
  SpSphericalTop::getParticleKineticEnergy(size_t ID) const {
    double KE = SpPoint::getParticleKineticEnergy(ID);

    const double I = getScalarMomentOfInertia(ID);
    if (!std::isinf(I))
      KE += 0.5 * I * Sim->dynamics->getRotData(ID).angularVelocity.nrm2();
    return KE;
  }

  void 
  SpSphericalTop::operator<<(const magnet::xml::Node& XML)
  {
    SpPoint::operator<<(XML);
    inertiaConstant = XML.getAttribute("InertiaConstant").as<double>() * Sim->units.unitInertia();
  }
}
