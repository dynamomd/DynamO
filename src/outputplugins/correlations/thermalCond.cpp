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

#include "thermalCond.hpp"
#include <boost/foreach.hpp>
#include "../../dynamics/include.hpp"
#include "../../dynamics/interactions/intEvent.hpp"
#include "../1partproperty/kenergy.hpp"
#include "../../base/is_ensemble.hpp"

COPThermalCon::COPThermalCon(const DYNAMO::SimData* tmp, const XMLNode& XML):
  COPCorrelator<CVector<> >(tmp,"ThermalConductivity", XML)
{
  
}

void 
COPThermalCon::initialise()
{
  COPCorrelator<CVector<> >::initialise();
  Sim->getOutputPlugin<COPKEnergy>();
  
  try {
    dynamic_cast<const DYNAMO::CENVE *>(Sim->Ensemble.get());
  }
  catch(std::exception)
    {
      D_throw() << "WARNING: This is only valid in the microcanonical"
	" ensemble!\nSee J.J. Erpenbeck, Phys. Rev. A 39, 4718 (1989) for more"
	"\n Essentially you need entropic data too for other ensembles";
    }

  dt = getdt();

  //Sum up the constant Del G.
  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    constDelG += part.getVelocity () * Sim->Dynamics.getParticleEnergy(part);

  I_cout() << "dt set to " << dt / Sim->Dynamics.units().unitTime();
}

Iflt 
COPThermalCon::rescaleFactor() 
{ 
  return Sim->Dynamics.units().unitk()
    / (Sim->Dynamics.units().unitTime() //This line should be 1 however we have scaled the correlator time as well
       * Sim->Dynamics.units().unitThermalCond() * 2.0 
       * count * pow(Sim->getOutputPlugin<COPKEnergy>()->getAvgkT(), 2)
       * Sim->Dynamics.units().simVolume());
}

inline void 
COPThermalCon::output(xmlw::XmlStream &XML)
{
  XML << xmlw::tag("EinsteinCorrelator")
      << xmlw::attr("name") << name
      << xmlw::attr("size") << accG2.size()
      << xmlw::attr("dt") << dt/Sim->Dynamics.units().unitTime()
      << xmlw::attr("LengthInMFT") << dt * accG2.size()
    / Sim->getOutputPlugin<COPMisc>()->getMFT()
      << xmlw::attr("simFactor") << rescaleFactor()
      << xmlw::attr("SampleCount") << count
      << xmlw::chardata();
  
  Iflt factor = rescaleFactor();
  
  for (unsigned int i = 0; i < accG2.size(); i++)
    {
      XML   << (1+i) * dt / Sim->Dynamics.units().unitTime()
	    << "\t ";
      
      for (int j=0;j<NDIM;j++)
	XML << accG2[i][j] * factor 
	    << "\t ";
      
      XML << "\n";
    }
  
  XML << xmlw::endtag("EinsteinCorrelator");
}

inline CVector<> 
COPThermalCon::impulseDelG(const C2ParticleData& PDat)
{
  return PDat.rij * PDat.particle1_.getDeltaeCalc();
}

void 
COPThermalCon::updateConstDelG(const C2ParticleData& PDat)
{
  Iflt p1E = Sim->Dynamics.getParticleEnergy(PDat.particle1_.getParticle());
  Iflt p2E = Sim->Dynamics.getParticleEnergy(PDat.particle2_.getParticle());
  
  constDelG += PDat.particle1_.getParticle().getVelocity() * p1E 
    + PDat.particle2_.getParticle().getVelocity() * p2E
    - PDat.particle1_.getOldVel() * (p1E - PDat.particle1_.getDeltaeCalc())
    - PDat.particle2_.getOldVel() * (p2E - PDat.particle2_.getDeltaeCalc());
}
