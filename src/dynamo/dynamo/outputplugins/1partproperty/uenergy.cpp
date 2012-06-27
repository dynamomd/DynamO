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

#include <dynamo/outputplugins/1partproperty/uenergy.hpp>
#include <dynamo/interactions/captures.hpp>
#include <dynamo/include.hpp>
#include <dynamo/interactions/intEvent.hpp>
#include <dynamo/simulation.hpp>
#include <boost/foreach.hpp>
#include <magnet/xmlwriter.hpp>
#include <cmath>

namespace dynamo {
  OPUEnergy::OPUEnergy(const dynamo::Simulation* tmp, const magnet::xml::Node&):
    OP1PP(tmp,"UEnergy", 250),
    intECurrent(0.0),
    intEsqAcc(0.0),
    intEAcc(0.0)
  {}

  void 
  OPUEnergy::changeSystem(OutputPlugin* Eplug)
  {
    std::swap(Sim, static_cast<OPUEnergy*>(Eplug)->Sim);
    std::swap(intECurrent, static_cast<OPUEnergy*>(Eplug)->intECurrent);
  }

  void
  OPUEnergy::initialise()
  {
    minE = maxE = intECurrent = Sim->calcInternalEnergy();
  }

  double 
  OPUEnergy::getAvgSqU() const
  { 
    return intEsqAcc / ( Sim->dSysTime 
			 * pow(Sim->units.unitEnergy(), 2)); 
  }

  double 
  OPUEnergy::getAvgU() const
  { 
    return intEAcc / ( Sim->dSysTime * Sim->units.unitEnergy()); 
  }

  void 
  OPUEnergy::A1ParticleChange(const ParticleEventData& PDat)
  {
    intECurrent += PDat.getDeltaU();
    minE = std::min(minE, intECurrent);
    maxE = std::max(maxE, intECurrent);
  }

  void 
  OPUEnergy::A2ParticleChange(const PairEventData& PDat)
  {
    intECurrent += PDat.particle1_.getDeltaU()
      + PDat.particle2_.getDeltaU();
    minE = std::min(minE, intECurrent);
    maxE = std::max(maxE, intECurrent);
  }

  void 
  OPUEnergy::stream(const double& dt)
  {
    intEAcc += intECurrent * dt;
    intEsqAcc += intECurrent * intECurrent * dt;
  }

  void
  OPUEnergy::output(magnet::xml::XmlStream &XML)
  {
    XML << magnet::xml::tag("CEnergy")
	<< magnet::xml::tag("InternalEnergy")
	<< magnet::xml::attr("Avg") << getAvgU()
	<< magnet::xml::attr("SquareAvg") << getAvgSqU()
	<< magnet::xml::attr("Current") << intECurrent
      / Sim->units.unitEnergy()
	<< magnet::xml::attr("Max") << maxE
      / Sim->units.unitEnergy()
	<< magnet::xml::attr("Min") << minE
      / Sim->units.unitEnergy()
	<< magnet::xml::endtag("InternalEnergy")
	<< magnet::xml::endtag("CEnergy");
  }

  void
  OPUEnergy::periodicOutput()
  {
    I_Pcout() << "U " << intECurrent / Sim->units.unitEnergy() << ", ";
  }
}
