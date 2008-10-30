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

#include "thermaldiffE.hpp"
#include <boost/foreach.hpp>
#include "../../dynamics/include.hpp"
#include "../../dynamics/interactions/intEvent.hpp"
#include "../1partproperty/kenergy.hpp"

COPThermalDiffusionE::COPThermalDiffusionE(const DYNAMO::SimData* tmp,
					 const XMLNode& XML):
  COPCorrelator<CVector<> >(tmp,"ThermalDiffusionE", XML),
  constDelGsp1(0.0),
  delGsp1(0.0),
  species1(0),
  sysMom(0.0),
  massFracSp1(1)
{
  operator<<(XML);
}

void 
COPThermalDiffusionE::operator<<(const XMLNode& XML)
{
  try 
    {
      try {
	species1 = Sim->Dynamics.getSpecies
	  (boost::lexical_cast<std::string>(XML.getAttribute("Species")))
	  .getID();
	
	COPCorrelator<CVector<> >::operator<<(XML);

      } catch (std::exception& nex)
	{
	  D_throw() << "The name of the Species must be specified"
		  << nex.what();
	}
    }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in COPMutualDiffusion";
    }
}

void 
COPThermalDiffusionE::initialise()
{
  try {
    dynamic_cast<const DYNAMO::CENVE *>(Sim->Ensemble.get());
  }
  catch(std::exception)
    {
      D_throw() << "WARNING: This is only valid in the microcanonical"
	" ensemble!\nSee J.J. Erpenbeck, Phys. Rev. A 39, 4718 (1989) for more"
	"\n Essentially you need entropic data too for other ensembles";
    }

  COPCorrelator<CVector<> >::initialise();
  Sim->getOutputPlugin<COPKEnergy>();
  
  accG2.resize(CorrelatorLength, CVector<>(0.0));
  Gsp1.resize(CorrelatorLength, CVector<>(0.0));
  dt = getdt();
  
  Iflt sysMass = 0.0;
  BOOST_FOREACH(const CSpecies& sp, Sim->Dynamics.getSpecies())
    sysMass += sp.getMass() * sp.getCount();

  //Sum up the constant Del G.
  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    {
      constDelG += part.getVelocity () * Sim->Dynamics.getParticleEnergy(part);
      sysMom += part.getVelocity() * Sim->Dynamics.getSpecies(part).getMass();
      
      if (Sim->Dynamics.getSpecies(part).getID() == species1)
	constDelGsp1 += part.getVelocity();
    }

  constDelGsp1 *= Sim->Dynamics.getSpecies()[species1].getMass();
  
  massFracSp1 = Sim->Dynamics.getSpecies()[species1].getCount() 
    * Sim->Dynamics.getSpecies()[species1].getMass() / sysMass;

  I_cout() << "dt set to " << dt / Sim->Dynamics.units().unitTime();
}

inline void 
COPThermalDiffusionE::output(xmlw::XmlStream &XML)
{
  XML << xmlw::tag("EinsteinCorrelator")
      << xmlw::attr("name") << name
      << xmlw::attr("size") << accG2.size()
      << xmlw::attr("dt") << dt/Sim->Dynamics.units().unitTime()
      << xmlw::attr("LengthInMFT") 
      << dt * accG2.size() / Sim->getOutputPlugin<COPMisc>()->getMFT()
      << xmlw::attr("simFactor") << rescaleFactor()
      << xmlw::attr("SampleCount") << count
      << xmlw::chardata();
  
  Iflt factor = rescaleFactor();
  
  for (size_t i = 0; i < accG2.size(); ++i)
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

Iflt 
COPThermalDiffusionE::rescaleFactor() 
{ 
  return 1.0
    / (Sim->Dynamics.units().unitTime() //This line should be 1 however we have scaled the correlator time as well
       * Sim->Dynamics.units().unitThermalDiffusion() * 2.0 
       * count * Sim->getOutputPlugin<COPKEnergy>()->getAvgkT()
       * Sim->Dynamics.units().simVolume());
}

void 
COPThermalDiffusionE::stream(const Iflt edt)
{      
  //Now test if we've gone over the step time
  if (currentdt + edt >= dt)
    {
      delG += constDelG * (dt - currentdt);
      delGsp1 += (constDelGsp1 - massFracSp1 * sysMom) * (dt - currentdt);
      newG();
      currentdt += edt - dt;
      
      while (currentdt >= dt)
	{
	  delG = constDelG * dt;
	  delGsp1 = (constDelGsp1 - massFracSp1 * sysMom) * dt;
	  currentdt -= dt;
	  newG();
	}

      //Now calculate the start of the new delG
      delG = constDelG * currentdt;
      delGsp1 = (constDelGsp1 - massFracSp1 * sysMom) * currentdt;
    }
  else
    {
      currentdt += edt;
      delG += constDelG * edt;
      delGsp1 += (constDelGsp1 - massFracSp1 * sysMom) * edt;
    }
}

void 
COPThermalDiffusionE::newG()
{
    //This ensures the list stays at accumilator size
  G.push_front (delG);
  Gsp1.push_front(delGsp1);

  if (notReady)
    {
      if (++currlen != CorrelatorLength)
	return;
      
      notReady = false;
    }
    
    accPass();
}

void 
COPThermalDiffusionE::accPass()
{
  ++count;
  CVector<> sum(0), sumsp1(0);

  for (size_t index = 0; index < CorrelatorLength; ++index)
    {
      sum += G[index];
      sumsp1 += Gsp1[index];
      accG2[index] += sum * sumsp1;
    }
}

inline CVector<> 
COPThermalDiffusionE::impulseDelG(const C2ParticleData& PDat)
{
  return PDat.rij * PDat.particle1_.getDeltaeCalc();
}

void 
COPThermalDiffusionE::updateConstDelG(const C2ParticleData& PDat)
{
  updateConstDelG(PDat.particle1_);
  updateConstDelG(PDat.particle2_);
}

void 
COPThermalDiffusionE::updateConstDelG(const C1ParticleData& PDat) 
{
  Iflt p1E = Sim->Dynamics.getParticleEnergy(PDat.getParticle());
  
  constDelG += PDat.getParticle().getVelocity() * p1E 
    - PDat.getOldVel() * (p1E - PDat.getDeltaeCalc());

  sysMom += PDat.getDeltaP();
  
  if (PDat.getSpecies().getID() == species1)
    constDelGsp1 += PDat.getDeltaP();
}
