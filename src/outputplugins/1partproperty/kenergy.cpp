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

#include "kenergy.hpp"
#include <boost/foreach.hpp>
#include <cmath>
#include "../../dynamics/interactions/captures.hpp"
#include "../../dynamics/include.hpp"
#include "../../dynamics/interactions/intEvent.hpp"
#include "../../base/is_simdata.hpp"
#include "../../extcode/xmlwriter.hpp"

OPKEnergy::OPKEnergy(const DYNAMO::SimData* tmp, const XMLNode&):
  OP1PP(tmp,"KEnergy", 250),
  InitialKE(0.0),
  KEacc(0.0),
  KEsqAcc(0.0),
  KECurrent(0.0)
{}

void
OPKEnergy::changeSystem(OutputPlugin* Eplug)
{
  std::swap(Sim, static_cast<OPKEnergy*>(Eplug)->Sim);

  std::swap(KECurrent,
	    static_cast<OPKEnergy*>(Eplug)->KECurrent);
}

void
OPKEnergy::temperatureRescale(const double& scale)
{
  KECurrent *= scale;
}

void
OPKEnergy::initialise()
{
  InitialKE = KECurrent = Sim->dynamics.getLiouvillean().getSystemKineticEnergy();
}

double
OPKEnergy::getAvgTheta() const
{
  return getAvgkT() / Sim->dynamics.units().unitEnergy();
}

double
OPKEnergy::getAvgkT() const
{
  return 2.0 * KEacc / (Sim->dSysTime * Sim->N * Sim->dynamics.getLiouvillean().getParticleDOF());
}

double
OPKEnergy::getAvgSqTheta() const
{
  return 2.0 * KEsqAcc / (Sim->dSysTime * Sim->N
			* Sim->dynamics.getLiouvillean().getParticleDOF()
			* Sim->dynamics.units().unitEnergy()
			* Sim->dynamics.units().unitEnergy());
}

void
OPKEnergy::A1ParticleChange(const ParticleEventData& PDat)
{
  KECurrent += PDat.getDeltaKE();
}

void
OPKEnergy::A2ParticleChange(const PairEventData& PDat)
{
  KECurrent += PDat.particle1_.getDeltaKE() + PDat.particle2_.getDeltaKE();
}

void
OPKEnergy::stream(const double& dt)
{
  KEacc += KECurrent * dt;
  KEsqAcc += KECurrent * KECurrent * dt;
}

void
OPKEnergy::output(xml::XmlStream &XML)
{
  double powerloss = (InitialKE - KECurrent) * Sim->dynamics.units().unitLength()
    * pow(Sim->dynamics.units().unitTime(),3)
    / (Sim->dynamics.units().unitMass() * Sim->dSysTime * Sim->dynamics.units().simVolume());

  XML << xml::tag("KEnergy")
      << xml::tag("T") << xml::attr("val") << getAvgTheta()
      << xml::attr("current")
      << (2.0 * Sim->dynamics.getLiouvillean().getSystemKineticEnergy()
	  / (Sim->dynamics.getLiouvillean().getParticleDOF() * Sim->N * Sim->dynamics.units().unitEnergy()))
      << xml::endtag("T")
      << xml::tag("T2") << xml::attr("val") << getAvgSqTheta()
      << xml::endtag("T2")

      << xml::tag("PowerLoss")
      << xml::attr("val") << powerloss
      << xml::endtag("PowerLoss")

      << xml::endtag("KEnergy");
}

void
OPKEnergy::periodicOutput()
{
  double powerloss = (InitialKE - KECurrent) * Sim->dynamics.units().unitLength()
    * pow(Sim->dynamics.units().unitTime(),3)
    / (Sim->dynamics.units().unitMass() * Sim->dSysTime * Sim->dynamics.units().simVolume());


  I_Pcout() << "T "
	    <<  2.0 * KECurrent / (Sim->N * Sim->dynamics.getLiouvillean().getParticleDOF()
				   * Sim->dynamics.units().unitEnergy())
	    << ", <T> " << getAvgTheta() << ", <PwrLoss> " << powerloss << ", ";
}

