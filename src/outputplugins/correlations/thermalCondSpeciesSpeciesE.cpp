/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "thermalCondSpeciesSpeciesE.hpp"
#include <boost/foreach.hpp>
#include "../../dynamics/include.hpp"
#include "../../dynamics/interactions/intEvent.hpp"
#include "../1partproperty/kenergy.hpp"
#include "../../base/is_ensemble.hpp"
#include "../0partproperty/misc.hpp"

OPThermalConductivitySpeciesSpeciesE::OPThermalConductivitySpeciesSpeciesE(const DYNAMO::SimData* tmp,
						 const XMLNode& XML):
  OutputPlugin(tmp,"ThermalConductivityE"),
  count(0),
  dt(0),
  currentdt(0.0),
  currlen(0),
  notReady(true),
  CorrelatorLength(100)
{
  operator<<(XML);
}

void 
OPThermalConductivitySpeciesSpeciesE::operator<<(const XMLNode& XML)
{
  try 
    {
      if (XML.isAttributeSet("Length"))
	CorrelatorLength = boost::lexical_cast<unsigned int>
	  (XML.getAttribute("Length"));
      
      if (XML.isAttributeSet("dt"))
	dt = Sim->dynamics.units().unitTime() * 
	  boost::lexical_cast<Iflt>(XML.getAttribute("dt"));
      
      if (XML.isAttributeSet("t"))
	dt = Sim->dynamics.units().unitTime() * 
	  boost::lexical_cast<Iflt>(XML.getAttribute("t"))/CorrelatorLength;
    }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in OPCorrelator";
    }
  
  }

void 
OPThermalConductivitySpeciesSpeciesE::initialise()
{
  size_t Nsp(Sim->dynamics.getSpecies().size());
  
  constDelG.resize(Nsp, Vector (0,0,0));
  
  delG.resize(Nsp, Vector (0,0,0));
  
  G.resize(Nsp, boost::circular_buffer<Vector  >(CorrelatorLength, Vector(0,0,0)));
  
  accG2.resize(Nsp * Nsp, std::vector<Vector  >(CorrelatorLength, Vector(0,0,0)));
  
  Sim->getOutputPlugin<OPMisc>();
  Sim->getOutputPlugin<OPKEnergy>();
  
  try {
    dynamic_cast<const DYNAMO::CENVE*>(Sim->Ensemble.get());
  }
  catch(std::exception)
    {
      D_throw() << "WARNING: This is only valid in the microcanonical"
	" ensemble!\nSee J.J. Erpenbeck, Phys. Rev. A 39, 4718 (1989) for more"
	"\n Essentially you need entropic data too for other ensembles";
    }

  if (dt == 0.0)
    {
      if (Sim->lastRunMFT != 0.0)
	dt = Sim->lastRunMFT * 50.0 / CorrelatorLength;
      else
	dt = 10.0 / (((Iflt) CorrelatorLength) 
		     * sqrt(Sim->dynamics.getLiouvillean().getkT()) * CorrelatorLength);
    }
  
  //Sum up the constant Del G.
  BOOST_FOREACH(const smrtPlugPtr<CSpecies>& spec, Sim->dynamics.getSpecies())
    BOOST_FOREACH(const size_t& id, *spec->getRange())
    {
      const CParticle& part(Sim->vParticleList[id]);
      constDelG[spec->getID()] += part.getVelocity() 
	* Sim->dynamics.getLiouvillean().getParticleKineticEnergy(part);
    }
  
  I_cout() << "dt set to " << dt / Sim->dynamics.units().unitTime();
}

Iflt 
OPThermalConductivitySpeciesSpeciesE::rescaleFactor() 
{ 
  return Sim->dynamics.units().unitk() 
    /(//This next line should be 1 however we have scaled the
      //correlator time as well
      Sim->dynamics.units().unitTime() 
      * Sim->dynamics.units().unitThermalCond() * 2.0 
      * count * pow(Sim->getOutputPlugin<OPKEnergy>()->getAvgkT(), 2)
      * Sim->dynamics.units().simVolume());
}

void 
OPThermalConductivitySpeciesSpeciesE::output(xmlw::XmlStream &XML)
{
  XML << xmlw::tag("EinsteinCorrelator")
      << xmlw::attr("name") << name
      << xmlw::attr("size") << accG2.size()
      << xmlw::attr("dt") << dt/Sim->dynamics.units().unitTime()
      << xmlw::attr("LengthInMFT") << dt * accG2.size()
    / Sim->getOutputPlugin<OPMisc>()->getMFT()
      << xmlw::attr("simFactor") << rescaleFactor()
      << xmlw::attr("SampleCount") << count;
  
  Iflt factor = rescaleFactor();

  size_t Nsp(Sim->dynamics.getSpecies().size());
  
  for (size_t id1(0); id1 < Nsp; ++id1)
    for (size_t id2(0); id2 < Nsp; ++id2)      
      {
	XML << xmlw::tag("Component") 
	    << xmlw::attr("Species1") << id1
	    << xmlw::attr("Species2") << id2
	    << xmlw::chardata();

	for (unsigned int i = 0; i < accG2.size(); i++)
	  {
	    XML << (1+i) * dt / Sim->dynamics.units().unitTime()
		<< "\t ";
	    
	    for (size_t j=0;j<NDIM;j++)
	      XML << accG2[id1 + Nsp * id2][i][j] * factor 
		  << "\t ";
	    
	    XML << "\n";
	  }

	XML << xmlw::endtag("Component");
      }
  
  XML << xmlw::endtag("EinsteinCorrelator");
}

void 
OPThermalConductivitySpeciesSpeciesE::stream(const Iflt& edt)
{
  size_t Nsp(Sim->dynamics.getSpecies().size());
  
  //Now test if we've gone over the step time
  if (currentdt + edt >= dt)
    {      
      for (size_t spec(0); spec < Nsp; ++spec)
	delG[spec] += constDelG[spec] * (dt - currentdt);
      
      newG();
      
      currentdt += edt - dt;
      
      while (currentdt >= dt)
	{
	  for (size_t spec(0); spec < Nsp; ++spec)
	    delG[spec] = constDelG[spec] * dt;
	  
	  currentdt -= dt;
	  newG();
	}
      
      //Now calculate the start of the new delG
      for (size_t spec(0); spec < Nsp; ++spec)
	delG[spec] = constDelG[spec] * currentdt;
    }
  else
    {
      currentdt += edt;

      for (size_t spec(0); spec < Nsp; ++spec)
	delG[spec] += constDelG[spec] * edt;
    }
}

void 
OPThermalConductivitySpeciesSpeciesE::eventUpdate(const CGlobEvent& iEvent, 
				     const CNParticleData& PDat) 
{
  stream(iEvent.getdt());
  //impulseDelG(PDat);
  updateConstDelG(PDat);
}

void 
OPThermalConductivitySpeciesSpeciesE::eventUpdate(const CLocalEvent& iEvent, 
				     const CNParticleData& PDat) 
{
  stream(iEvent.getdt());
  //impulseDelG(PDat);
  updateConstDelG(PDat);
}

void 
OPThermalConductivitySpeciesSpeciesE::eventUpdate(const CSystem&, 
				     const CNParticleData& PDat,
				     const Iflt& edt) 
{ 
  stream(edt);
  //impulseDelG(PDat);
  updateConstDelG(PDat);
}
  
void 
OPThermalConductivitySpeciesSpeciesE::eventUpdate(const CIntEvent& iEvent,
				     const C2ParticleData& PDat)
{
  stream(iEvent.getdt());
  //impulseDelG(PDat);
  updateConstDelG(PDat);
}

void
OPThermalConductivitySpeciesSpeciesE::impulseDelG(const CNParticleData& ndat) 
{ 
  /*  Vector  acc(0);
  
  BOOST_FOREACH(const C2ParticleData& dat, ndat.L2partChanges)
    acc += impulseDelG(dat);
  
    return acc;*/
}

void
OPThermalConductivitySpeciesSpeciesE::impulseDelG(const C2ParticleData& PDat)
{
  //return PDat.rij * PDat.particle1_.getDeltaeCalc();
}

void 
OPThermalConductivitySpeciesSpeciesE::newG()
{
  //This ensures the list stays at accumilator size
  
  size_t Nsp(Sim->dynamics.getSpecies().size());
  
  for (size_t id(0); id < Nsp; ++id)
    G[id].push_front(delG[id]);
  
  if (notReady)
    {
      if (++currlen != CorrelatorLength)
	return;
      
      notReady = false;
    }
  
  accPass();
}

void 
OPThermalConductivitySpeciesSpeciesE::accPass()
{
  ++count;

  size_t Nsp(Sim->dynamics.getSpecies().size());
  
  for (size_t id1(0); id1 < Nsp; ++id1)
    for (size_t id2(0); id2 < Nsp; ++id2)
      {
	Vector  sum1(0,0,0), sum2(0,0,0);
	
	for (size_t i = 0; i < CorrelatorLength; ++i)
	  {
	    sum1 += G[id1][i];
	    sum2 += G[id2][i];

	    Vector  tmp (sum1);

	    for (size_t j(0); j < NDIM; ++j)
	      tmp[j] *= sum2[j];
		
	    accG2[id1+Nsp*id2][i] += tmp;
	  }
      }
}

void 
OPThermalConductivitySpeciesSpeciesE::updateConstDelG(const CNParticleData& ndat)
{
  BOOST_FOREACH(const C1ParticleData& dat, ndat.L1partChanges)
    updateConstDelG(dat);
  
  BOOST_FOREACH(const C2ParticleData& dat, ndat.L2partChanges)
    updateConstDelG(dat);
}

void 
OPThermalConductivitySpeciesSpeciesE::updateConstDelG(const C1ParticleData& PDat)
{
  Iflt p1E = Sim->dynamics.getLiouvillean().getParticleKineticEnergy(PDat.getParticle());
  
  constDelG[PDat.getSpecies().getID()] += PDat.getParticle().getVelocity() * p1E 
    - PDat.getOldVel() * (p1E - PDat.getDeltaKE());
}

void 
OPThermalConductivitySpeciesSpeciesE::updateConstDelG(const C2ParticleData& PDat)
{
  Iflt p1E = Sim->dynamics.getLiouvillean().getParticleKineticEnergy(PDat.particle1_.getParticle());
  Iflt p2E = Sim->dynamics.getLiouvillean().getParticleKineticEnergy(PDat.particle2_.getParticle());
  
  constDelG[Sim->dynamics.getSpecies(PDat.particle1_.getParticle()).getID()]
    += PDat.particle1_.getParticle().getVelocity() * p1E
    - PDat.particle1_.getOldVel() * (p1E - PDat.particle1_.getDeltaKE());
  
  constDelG[Sim->dynamics.getSpecies(PDat.particle2_.getParticle()).getID()]
    += PDat.particle2_.getParticle().getVelocity() * p2E
    - PDat.particle2_.getOldVel() * (p2E - PDat.particle2_.getDeltaKE());
}
