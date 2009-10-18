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

#include "kenergy.hpp"
#include <boost/foreach.hpp>
#include <cmath>
#include "../../dynamics/interactions/captures.hpp"
#include "../../dynamics/include.hpp"
#include "../../dynamics/interactions/intEvent.hpp"
#include "../../base/is_simdata.hpp"

COPKEnergy::COPKEnergy(const DYNAMO::SimData* tmp, const XMLNode&):
  COP1PP(tmp,"KEnergy", 250),
  InitialKE(0.0),
  KEacc(0.0),
  KEsqAcc(0.0),
  KECurrent(0.0)
{}

void 
COPKEnergy::changeSystem(COutputPlugin* Eplug)
{
  std::swap(Sim, static_cast<COPKEnergy*>(Eplug)->Sim);
  
  std::swap(KECurrent, 
	    static_cast<COPKEnergy*>(Eplug)->KECurrent);
}

void 
COPKEnergy::temperatureRescale(const Iflt& scale)
{
  KECurrent *= scale * scale;
}

void
COPKEnergy::initialise()
{  
  InitialKE = KECurrent = Sim->Dynamics.getLiouvillean().getSystemKineticEnergy();
}

Iflt 
COPKEnergy::getAvgTheta() const
{
  return getAvgkT() / Sim->Dynamics.units().unitEnergy();
}

Iflt 
COPKEnergy::getAvgkT() const
{
  return 2.0 * KEacc / (Sim->dSysTime * Sim->lN * Sim->Dynamics.getLiouvillean().getParticleDOF());
}

Iflt 
COPKEnergy::getAvgSqTheta() const
{
  return 2.0 * KEsqAcc / (Sim->dSysTime * Sim->lN 
			* Sim->Dynamics.getLiouvillean().getParticleDOF()
			* Sim->Dynamics.units().unitEnergy()
			* Sim->Dynamics.units().unitEnergy());
}

void 
COPKEnergy::A1ParticleChange(const C1ParticleData& PDat)
{
  KECurrent += PDat.getDeltaKE();
}

void
COPKEnergy::A2ParticleChange(const C2ParticleData& PDat)
{
  KECurrent += PDat.particle1_.getDeltaKE() + PDat.particle2_.getDeltaKE();
}

void 
COPKEnergy::stream(const Iflt& dt)
{
  KEacc += KECurrent * dt;
  KEsqAcc += KECurrent * KECurrent * dt;
}

void
COPKEnergy::output(xmlw::XmlStream &XML)
{
  Iflt powerloss = (InitialKE - KECurrent) * Sim->Dynamics.units().unitLength()
    * pow(Sim->Dynamics.units().unitTime(),3) 
    / (Sim->Dynamics.units().unitMass() * Sim->dSysTime * Sim->Dynamics.units().simVolume());

  XML << xmlw::tag("KEnergy")
      << xmlw::tag("T") << xmlw::attr("val") << getAvgTheta()
      << xmlw::attr("current") 
      << (2.0 * Sim->Dynamics.getLiouvillean().getSystemKineticEnergy() 
	  / (static_cast<Iflt>(NDIM) * Sim->lN * Sim->Dynamics.units().unitEnergy()))
      << xmlw::endtag("T")
      << xmlw::tag("T2") << xmlw::attr("val") << getAvgSqTheta()
      << xmlw::endtag("T2")
    
      << xmlw::tag("PowerLoss")
      << xmlw::attr("val") << powerloss
      << xmlw::endtag("PowerLoss")

      << xmlw::endtag("KEnergy");
}

void
COPKEnergy::periodicOutput()
{
  Iflt powerloss = (InitialKE - KECurrent) * Sim->Dynamics.units().unitLength()
    * pow(Sim->Dynamics.units().unitTime(),3) 
    / (Sim->Dynamics.units().unitMass() * Sim->dSysTime * Sim->Dynamics.units().simVolume());


  I_Pcout() << "T " 
	    <<  2.0 * KECurrent / (Sim->lN * Sim->Dynamics.getLiouvillean().getParticleDOF() 
				   * Sim->Dynamics.units().unitEnergy())
	    << ", <T> " << getAvgTheta() << ", <PwrLoss> " << powerloss << ", ";
}

