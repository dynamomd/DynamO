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

#include "thermaldiff.hpp"
#include <boost/foreach.hpp>
#include "../../dynamics/include.hpp"
#include "../../dynamics/interactions/intEvent.hpp"

COPThermalDiffusion::COPThermalDiffusion(const DYNAMO::SimData* tmp,
					 const XMLNode& XML):
  COPCorrelator<CVector<> >(tmp,"ThermalDiffusion", XML),
  constDelGsp1(0.0),
  delGsp1(0.0),
  species1(NULL),
  sysMom(0.0),
  massFracSp1(1)
{  
  try 
    {
      try {
	species1 = &Sim->Dynamics.getSpecies(boost::lexical_cast<std::string>(XML.getAttribute("Species")));
	
	COPCorrelator<CVector<> >::operator<<(XML);
      } catch (std::exception& nex)
	{
	  D_throw() << "Failed to find the species for the mutual diffusion\n"
		  << nex.what();
	}
    }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in COPMutualDiffusion";
    }
}


void 
COPThermalDiffusion::initialise()
{
  COPCorrelator<CVector<> >::initialise();
  accG2.resize(CorrelatorLength, CVector<>(0.0));
  dt = getdt();
  
  Iflt sysMass = 0.0;
  //Sum up the constant Del G.
  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    {
      constDelG += part.getVelocity () * Sim->Dynamics.getParticleEnergy(part);
      
      sysMom += part.getVelocity() * Sim->Dynamics.getSpecies(part).getMass();
      sysMass += Sim->Dynamics.getSpecies(part).getMass();
      
      if (species1->isSpecies(part))
	constDelGsp1 += part.getVelocity();      
    }

  constDelGsp1 *= species1->getMass();
  
  massFracSp1 = species1->getCount() * species1->getMass() / sysMass; 

  I_cout() << "WARNING: This is only valid in the microcanonical ensemble!";
}

inline void 
COPThermalDiffusion::output(xmlw::XmlStream &XML)
{
  XML << xmlw::tag("EinsteinCorrelator")
      << xmlw::attr("name") << name
      << xmlw::attr("size") << accG2.size()
      << xmlw::attr("dt") << dt/Sim->Dynamics.units().unitTime()
      << xmlw::attr("LengthInMFT") << dt * accG2.size() / ptrMisc->getMFT()
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

Iflt 
COPThermalDiffusion::rescaleFactor() 
{ 
  return 1.0
    / (Sim->Dynamics.units().unitTime() //This line should be 1 however we have scaled the correlator time as well
       * Sim->Dynamics.units().unitThermalDiffusion() * 2.0 
       * count * ptrEnergy->getAvgkT()
       * Sim->Dynamics.units().simVolume());
}

void 
COPThermalDiffusion::stream(const Iflt edt)
{    
  //Move the time forward
  currentdt += edt;
  
  //Now test if we've gone over the step time
  if (currentdt >= dt)
    {
      currentdt -= dt;
      delG += constDelG * (edt - currentdt);
      delGsp1 += (constDelGsp1 - massFracSp1 * sysMom) * (edt - currentdt);
      newG();
      
      //Now calculate the start of the new delG
      delG = constDelG * currentdt;
      delGsp1 = (constDelGsp1 - massFracSp1 * sysMom) * currentdt;
    }
  else
    {
      delG += constDelG * edt;
      delGsp1 += (constDelGsp1 - massFracSp1 * sysMom) * edt;
    }
}

void 
COPThermalDiffusion::newG()
{
  //This ensures the list stays at accumilator size
  if (G.size () == accG2.size ())
    {
      G.pop_back ();
      G.push_front (delG);
      
      Gsp1.pop_back();
      Gsp1.push_front(delGsp1);

      accPass ();
    }
  else
    {
      G.push_front (delG);
      Gsp1.push_front (delGsp1);

      if (G.size () == accG2.size ())
	accPass();
    }
}

void 
COPThermalDiffusion::accPass()
{
  count++;
  CVector<> sum(0), sumsp1(0);
  
  std::list<CVector<> >::iterator iptr = G.begin();
  std::list<CVector<> >::iterator iptrsp1 = Gsp1.begin();

  for (unsigned int index = 0; index < accG2.size(); index++)
    {
      sum += *iptr;
      sumsp1 += *iptrsp1;
      accG2[index] += sum * sumsp1;
      iptr++;
      iptrsp1++;
    }
}

inline CVector<> 
COPThermalDiffusion::impulseDelG(const C2ParticleData& PDat)
{
  return PDat.rij * PDat.particle1_.getDeltaeCalc();
}

void 
COPThermalDiffusion::updateConstDelG(const C2ParticleData& PDat)
{
  updateConstDelG(PDat.particle1_);
  updateConstDelG(PDat.particle2_);
}

void 
COPThermalDiffusion::updateConstDelG(const C1ParticleData& PDat) 
{
  Iflt p1E = Sim->Dynamics.getParticleEnergy(PDat.getParticle());
  
  constDelG += PDat.getParticle().getVelocity() * p1E 
    - PDat.getOldVel() * (p1E - PDat.getDeltaeCalc());

  sysMom += PDat.getDeltaP();
  
  if (species1->isSpecies(PDat.getParticle()))
    constDelGsp1 += PDat.getDeltaP();
}
