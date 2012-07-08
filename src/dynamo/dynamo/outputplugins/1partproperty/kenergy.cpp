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

#include <dynamo/outputplugins/1partproperty/kenergy.hpp>
#include <dynamo/interactions/captures.hpp>
#include <dynamo/include.hpp>
#include <dynamo/interactions/intEvent.hpp>
#include <dynamo/simulation.hpp>
#include <boost/foreach.hpp>
#include <magnet/xmlwriter.hpp>
#include <cmath>

namespace dynamo {
  OPKEnergy::OPKEnergy(const dynamo::Simulation* tmp, const magnet::xml::Node&):
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
    InitialKE = KECurrent = Sim->dynamics->getSystemKineticEnergy();
  }

  double
  OPKEnergy::getAvgTheta() const
  {
    return getAvgkT() / Sim->units.unitEnergy();
  }

  double
  OPKEnergy::getAvgkT() const
  {
    return 2.0 * KEacc / (Sim->systemTime * Sim->N * Sim->dynamics->getParticleDOF());
  }

  double
  OPKEnergy::getAvgSqTheta() const
  {
    return 2.0 * KEsqAcc / (Sim->systemTime * Sim->N
			    * Sim->dynamics->getParticleDOF()
			    * Sim->units.unitEnergy()
			    * Sim->units.unitEnergy());
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
  OPKEnergy::output(magnet::xml::XmlStream &XML)
  {
    double powerloss = (InitialKE - KECurrent) * Sim->units.unitLength()
      * pow(Sim->units.unitTime(),3)
      / (Sim->units.unitMass() * Sim->systemTime * Sim->getSimVolume());

    XML << magnet::xml::tag("KEnergy")
	<< magnet::xml::tag("T") << magnet::xml::attr("val") << getAvgTheta()
	<< magnet::xml::attr("current")
	<< (2.0 * Sim->dynamics->getSystemKineticEnergy()
	    / (Sim->dynamics->getParticleDOF() * Sim->N * Sim->units.unitEnergy()))
	<< magnet::xml::endtag("T")
	<< magnet::xml::tag("T2") << magnet::xml::attr("val") << getAvgSqTheta()
	<< magnet::xml::endtag("T2")

	<< magnet::xml::tag("PowerLoss")
	<< magnet::xml::attr("val") << powerloss
	<< magnet::xml::endtag("PowerLoss")

	<< magnet::xml::endtag("KEnergy");
  }

  void
  OPKEnergy::periodicOutput()
  {
    double powerloss = (InitialKE - KECurrent) * Sim->units.unitLength()
      * pow(Sim->units.unitTime(),3)
      / (Sim->units.unitMass() * Sim->systemTime * Sim->getSimVolume());


    I_Pcout() << "T "
	      <<  2.0 * KECurrent / (Sim->N * Sim->dynamics->getParticleDOF()
				     * Sim->units.unitEnergy())
	      << ", <T> " << getAvgTheta() << ", <PwrLoss> " << powerloss << ", ";
  }
}
