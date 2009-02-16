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

#include "kenergy.hpp"
#include <boost/foreach.hpp>
#include <cmath>
#include "../../dynamics/interactions/captures.hpp"
#include "../../dynamics/include.hpp"
#include "../../dynamics/interactions/intEvent.hpp"
#include "../../base/is_simdata.hpp"

COPKEnergy::COPKEnergy(const DYNAMO::SimData* tmp, const XMLNode&):
  COP1PP(tmp,"KEnergy", 250),
  sumPowerLoss(0.0),
  kTacc(0.0),
  kTsqAcc(0.0),
  kTCurrent(0.0)
{}

void 
COPKEnergy::changeSystem(COutputPlugin* Eplug)
{
  std::swap(Sim, static_cast<COPKEnergy*>(Eplug)->Sim);

  for (int iDim = 0; iDim < NDIM; ++iDim)
    std::swap(kTCurrent[iDim], 
	      static_cast<COPKEnergy*>(Eplug)->kTCurrent[iDim]);
}

void 
COPKEnergy::temperatureRescale(const Iflt& scale)
{
  for (int iDim = 0; iDim < NDIM; ++iDim)
    kTCurrent[iDim] *= scale * scale;
}

void
COPKEnergy::initialise()
{  
  kTCurrent = Sim->Dynamics.getVeckT();
}

Iflt 
COPKEnergy::getAvgTheta() const
{
  return getAvgkT() / Sim->Dynamics.units().unitEnergy();
}

Iflt 
COPKEnergy::getAvgkT() const
{
  Iflt sum = 0.0;
  
  for (int iDim = 0; iDim < NDIM; iDim++)
    sum += kTacc[iDim];
  
  sum /= ((Iflt) NDIM) * Sim->dSysTime;
  
  return sum;
}

Iflt 
COPKEnergy::getAvgTheta(long i) const
{ 
  return kTacc[i] / (Sim->dSysTime * Sim->Dynamics.units().unitEnergy()); 
}


Iflt 
COPKEnergy::getAvgSqTheta() const
{
  Iflt sum = 0.0;
  
  for (int iDim = 0; iDim < NDIM; iDim++)
    sum += kTsqAcc[iDim];
  
  sum /= ((Iflt) NDIM) * Sim->dSysTime 
    * pow(Sim->Dynamics.units().unitEnergy(),2);
  
  return sum;
}

Iflt 
COPKEnergy::getAvgSqTheta(long i) const
{ 
  return kTsqAcc[i]/( Sim->dSysTime 
		      * pow(Sim->Dynamics.units().unitEnergy(),2)); 
}

void 
COPKEnergy::A1ParticleChange(const C1ParticleData& PDat)
{
  sumPowerLoss += PDat.getDeltae();
  //Update the accumilators
  kTCurrent += PDat.deltake * 2.0 / Sim->lN;    
}

void 
COPKEnergy::A2ParticleChange(const C2ParticleData& PDat)
{
  sumPowerLoss += PDat.getDeltae();
  //Update the accumilators
  kTCurrent += PDat.deltake * 2.0 / Sim->lN;  
}

void 
COPKEnergy::stream(const Iflt& dt)
{
  kTacc += kTCurrent * dt;
  kTsqAcc += kTCurrent * kTCurrent * dt;
}

void
COPKEnergy::output(xmlw::XmlStream &XML)
{
  Iflt unitPwrloss = Sim->Dynamics.units().unitLength()
    * pow(Sim->Dynamics.units().unitTime(),3) / Sim->Dynamics.units().unitMass();

  XML << xmlw::tag("KEnergy")
      << xmlw::tag("T") << xmlw::attr("val") << getAvgTheta()
      << xmlw::attr("current") << (2.0 * Sim->Dynamics.getKineticEnergy() / (static_cast<Iflt>(NDIM) * Sim->lN * Sim->Dynamics.units().unitEnergy()))
      << xmlw::endtag("T")
      << xmlw::tag("T2") << xmlw::attr("val") << getAvgSqTheta()
      << xmlw::endtag("T2")
    
      << xmlw::tag("Tx") 
      << xmlw::attr("val") << getAvgTheta(0) 
      << xmlw::endtag("Tx")
      << xmlw::tag("Ty") 
      << xmlw::attr("val") << getAvgTheta(1) 
      << xmlw::endtag("Ty")
      << xmlw::tag("Tz") 
      << xmlw::attr("val") << getAvgTheta(2) 
      << xmlw::endtag("Tz")

      << xmlw::tag("T2x") 
      << xmlw::attr("val") << getAvgSqTheta(0) 
      << xmlw::endtag("T2x")
      << xmlw::tag("T2y") 
      << xmlw::attr("val") << getAvgSqTheta(1) 
      << xmlw::endtag("T2y")
      << xmlw::tag("T2z") 
      << xmlw::attr("val") << getAvgSqTheta(2) 
      << xmlw::endtag("T2z")

      << xmlw::tag("PowerLoss")
      << xmlw::attr("val") << (sumPowerLoss*unitPwrloss/Sim->dSysTime)
    /Sim->Dynamics.units().simVolume()
      << xmlw::endtag("PowerLoss")

      << xmlw::endtag("KEnergy");
}

void
COPKEnergy::periodicOutput()
{
  Iflt unitPwrloss = Sim->Dynamics.units().unitLength()
    * pow(Sim->Dynamics.units().unitTime(),3) / Sim->Dynamics.units().unitMass();

  Iflt kT = 0.0;

  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    kT += kTCurrent[iDim];

  kT /= NDIM;

  I_Pcout() << "T " <<  kT / Sim->Dynamics.units().unitEnergy() << ", <T> " 
	    << getAvgTheta() << ", <PwrLoss> " 
	    << sumPowerLoss * unitPwrloss
    /(Sim->Dynamics.units().simVolume()*Sim->dSysTime)
	    << ", ";
}

