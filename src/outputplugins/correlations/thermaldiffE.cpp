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
  COutputPlugin(tmp,"ThermalDiffusionE", 60),
  G(100),
  count(0),
  dt(0),
  currentdt(0.0),
  constDelG(0.0), 
  delG(0.0),
  currlen(0),
  notReady(true),
  CorrelatorLength(100),
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
	species1name = boost::lexical_cast<std::string>
	  (XML.getAttribute("Species"));
	  
      } catch (std::exception& nex)
	{
	  D_throw() << "The name of the Species must be specified"
		    << nex.what();
	}
      
      if (XML.isAttributeSet("Length"))
	CorrelatorLength = boost::lexical_cast<unsigned int>
	  (XML.getAttribute("Length"));
      
      if (XML.isAttributeSet("dt"))
	dt = Sim->Dynamics.units().unitTime() * 
	  boost::lexical_cast<Iflt>(XML.getAttribute("dt"));
      
      if (XML.isAttributeSet("t"))
	dt = Sim->Dynamics.units().unitTime() * 
	  boost::lexical_cast<Iflt>(XML.getAttribute("t"))/CorrelatorLength;
    }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in COPMutualDiffusion";
    }
}

void 
COPThermalDiffusionE::initialise()
{
  species1 = Sim->Dynamics.getSpecies(species1name).getID();

  try {
    dynamic_cast<const DYNAMO::CENVE *>(Sim->Ensemble.get());
  }
  catch(std::exception)
    {
      D_throw() << "WARNING: This is only valid in the microcanonical"
	" ensemble!\nSee J.J. Erpenbeck, Phys. Rev. A 39, 4718 (1989) for more"
	"\n Essentially you need entropic data too for other ensembles";
    }

  G.resize(CorrelatorLength, CVector<>(0.0));
  accG2.resize(CorrelatorLength, CVector<>(0.0));
  Sim->getOutputPlugin<COPMisc>();
  Sim->getOutputPlugin<COPKEnergy>();
  
  accG2.resize(CorrelatorLength, CVector<>(0.0));
  Gsp1.resize(CorrelatorLength, CVector<>(0.0));

  if (dt == 0.0)
    {
      if (Sim->lastRunMFT != 0.0)
	dt = Sim->lastRunMFT * 50.0 / CorrelatorLength;
      else
	dt = 10.0 / (((Iflt) CorrelatorLength) 
		     * sqrt(Sim->Dynamics.Liouvillean().getkT()) * CorrelatorLength);
    }
  
  Iflt sysMass = 0.0;
  BOOST_FOREACH(const CSpecies& sp, Sim->Dynamics.getSpecies())
    sysMass += sp.getMass() * sp.getCount();

  //Sum up the constant Del G.
  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    {
      constDelG += part.getVelocity () * Sim->Dynamics.Liouvillean().getParticleKineticEnergy(part);
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
    / (Sim->Dynamics.units().unitTime() 
       // /\ /\ This line should be 1 however we have scaled the
       //correlator time as well
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
      
      CVector<> tmp(sum);
      
      for (size_t i(0); i < NDIM; ++i)
	tmp[i] *= sumsp1[i];

      accG2[index] += tmp;
    }
}

inline CVector<> 
COPThermalDiffusionE::impulseDelG(const C2ParticleData& PDat)
{
  return PDat.rij * PDat.particle1_.getDeltaKE();
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
  Iflt p1E = Sim->Dynamics.Liouvillean().getParticleKineticEnergy(PDat.getParticle());
  
  constDelG += PDat.getParticle().getVelocity() * p1E 
    - PDat.getOldVel() * (p1E - PDat.getDeltaKE());

  sysMom += PDat.getDeltaP();
  
  if (PDat.getSpecies().getID() == species1)
    constDelGsp1 += PDat.getDeltaP();
}

void 
COPThermalDiffusionE::eventUpdate(const CGlobEvent& iEvent, 
				  const CNParticleData& PDat) 
{
  stream(iEvent.getdt());
  delG += impulseDelG(PDat);
  updateConstDelG(PDat);
}

void 
COPThermalDiffusionE::eventUpdate(const CLocalEvent& iEvent, 
				  const CNParticleData& PDat) 
{
  stream(iEvent.getdt());
  delG += impulseDelG(PDat);
  updateConstDelG(PDat);
}

void 
COPThermalDiffusionE::eventUpdate(const CSystem&, 
				  const CNParticleData& PDat, 
				  const Iflt& edt) 
{ 
  stream(edt);
  delG += impulseDelG(PDat);
  updateConstDelG(PDat);
}
  
void 
COPThermalDiffusionE::eventUpdate(const CIntEvent& iEvent, 
				  const C2ParticleData& PDat)
{
  stream(iEvent.getdt());
  delG += impulseDelG(PDat);
  updateConstDelG(PDat);
}

CVector<> 
COPThermalDiffusionE::impulseDelG(const CNParticleData& ndat) 
{ 
  CVector<> acc(0);
    
  BOOST_FOREACH(const C2ParticleData& dat, ndat.L2partChanges)
    acc += impulseDelG(dat);
  
  return acc;
}

void 
COPThermalDiffusionE::updateConstDelG(const CNParticleData& ndat)
{
  BOOST_FOREACH(const C1ParticleData& dat, ndat.L1partChanges)
    updateConstDelG(dat);
  
  BOOST_FOREACH(const C2ParticleData& dat, ndat.L2partChanges)
    updateConstDelG(dat);
}
