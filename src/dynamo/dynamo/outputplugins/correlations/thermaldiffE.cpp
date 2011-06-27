/*  dynamo:- Event driven molecular dynamics simulator 
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

#include "thermaldiffE.hpp"
#include "../../dynamics/include.hpp"
#include "../../dynamics/interactions/intEvent.hpp"
#include "../1partproperty/kenergy.hpp"
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>

OPThermalDiffusionE::OPThermalDiffusionE(const dynamo::SimData* tmp,
					 const magnet::xml::Node& XML):
  OutputPlugin(tmp,"ThermalDiffusionE", 60),
  G(100),
  count(0),
  dt(0),
  currentdt(0.0),
  constDelG(0,0,0), 
  delG(0,0,0),
  currlen(0),
  notReady(true),
  CorrelatorLength(100),
  constDelGsp1(0,0,0),
  delGsp1(0,0,0),
  species1(0),
  sysMom(0,0,0),
  massFracSp1(1)
{
  operator<<(XML);
}

void 
OPThermalDiffusionE::operator<<(const magnet::xml::Node& XML)
{
  try 
    {
      try {
	species1name = std::string(XML.getAttribute("Species"));
	
      } catch (std::exception& nex)
	{
	  M_throw() << "The name of the Species must be specified"
		    << nex.what();
	}

      CorrelatorLength = XML.getAttribute("Length").as<size_t>(100);
      
      if (XML.getAttribute("dt").valid())
	dt = Sim->dynamics.units().unitTime() * 
	  XML.getAttribute("dt").as<double>();
      
      if (XML.getAttribute("t").valid())
	dt = Sim->dynamics.units().unitTime() * 
	  XML.getAttribute("t").as<double>() / CorrelatorLength;
    }
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in OPMutualDiffusion";
    }
}

void 
OPThermalDiffusionE::initialise()
{
  species1 = Sim->dynamics.getSpecies(species1name).getID();

  if (dynamic_cast<const dynamo::EnsembleNVE* >(Sim->ensemble.get()) == NULL)
    M_throw() << "WARNING: This is only valid in the microcanonical"
      " ensemble!\nSee J.J. Erpenbeck, Phys. Rev. A 39, 4718 (1989) for more"
      "\n Essentially you need entropic data too for other ensembles";

  G.resize(CorrelatorLength, Vector (0,0,0));
  accG2.resize(CorrelatorLength, Vector (0,0,0));
  Sim->getOutputPlugin<OPMisc>();
  Sim->getOutputPlugin<OPKEnergy>();
  
  accG2.resize(CorrelatorLength, Vector (0,0,0));
  Gsp1.resize(CorrelatorLength, Vector (0,0,0));

  if (dt == 0.0)
    {
      if (Sim->lastRunMFT != 0.0)
	dt = Sim->lastRunMFT * 50.0 / CorrelatorLength;
      else
	dt = 10.0 / (((double) CorrelatorLength) 
		     * sqrt(Sim->dynamics.getLiouvillean().getkT()) * CorrelatorLength);
    }
  
  double sysMass = 0;
  BOOST_FOREACH(const magnet::ClonePtr<Species>& sp, Sim->dynamics.getSpecies())
    BOOST_FOREACH(const size_t ID, *(sp->getRange()))
    sysMass += sp->getMass(ID);

  //Sum up the constant Del G and the mass fraction of the species
  double speciesMass = 0;
  BOOST_FOREACH(const Particle& part, Sim->particleList)
    {
      double mass = Sim->dynamics.getSpecies(part).getMass(part.getID());

      constDelG += part.getVelocity() * mass
	* Sim->dynamics.getLiouvillean().getParticleKineticEnergy(part);
      sysMom += part.getVelocity() * mass;
      
      if (Sim->dynamics.getSpecies(part).getID() == species1)
	{
	  constDelGsp1 += part.getVelocity();
	  speciesMass += mass;
	}
    }

  massFracSp1 = speciesMass / sysMass;
  dout << "dt set to " << dt / Sim->dynamics.units().unitTime() << std::endl;
}

inline void 
OPThermalDiffusionE::output(magnet::xml::XmlStream &XML)
{
  XML << magnet::xml::tag("EinsteinCorrelator")
      << magnet::xml::attr("name") << name
      << magnet::xml::attr("size") << accG2.size()
      << magnet::xml::attr("dt") << dt/Sim->dynamics.units().unitTime()
      << magnet::xml::attr("LengthInMFT") 
      << dt * accG2.size() / Sim->getOutputPlugin<OPMisc>()->getMFT()
      << magnet::xml::attr("simFactor") << rescaleFactor()
      << magnet::xml::attr("SampleCount") << count
      << magnet::xml::chardata();
  
  double factor = rescaleFactor();
  
  for (size_t i = 0; i < accG2.size(); ++i)
    {
      XML   << (1+i) * dt / Sim->dynamics.units().unitTime()
	    << "\t ";
      
      for (size_t j=0;j<NDIM;j++)
	XML << accG2[i][j] * factor 
	    << "\t ";
      
      XML << "\n";
    }
  
  XML << magnet::xml::endtag("EinsteinCorrelator");
}

double 
OPThermalDiffusionE::rescaleFactor() 
{ 
  return 1.0
    / (Sim->dynamics.units().unitTime() 
       // /\ /\ This line should be 1 however we have scaled the
       //correlator time as well
       * Sim->dynamics.units().unitThermalDiffusion() * 2.0 
       * count * Sim->getOutputPlugin<OPKEnergy>()->getAvgkT()
       * Sim->dynamics.getSimVolume());
}

void 
OPThermalDiffusionE::stream(const double edt)
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
OPThermalDiffusionE::newG()
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
OPThermalDiffusionE::accPass()
{
  ++count;
  Vector  sum(0,0,0), sumsp1(0,0,0);

  for (size_t index = 0; index < CorrelatorLength; ++index)
    {
      sum += G[index];
      sumsp1 += Gsp1[index];
      
      Vector  tmp(sum);
      
      for (size_t i(0); i < NDIM; ++i)
	tmp[i] *= sumsp1[i];

      accG2[index] += tmp;
    }
}

inline Vector  
OPThermalDiffusionE::impulseDelG(const PairEventData& PDat)
{
  return PDat.rij * PDat.particle1_.getDeltaKE();
}

void 
OPThermalDiffusionE::updateConstDelG(const PairEventData& PDat)
{
  updateConstDelG(PDat.particle1_);
  updateConstDelG(PDat.particle2_);
}

void 
OPThermalDiffusionE::updateConstDelG(const ParticleEventData& PDat) 
{
  double p1E = Sim->dynamics.getLiouvillean().getParticleKineticEnergy(PDat.getParticle());
  
  constDelG += PDat.getParticle().getVelocity() * p1E 
    - PDat.getOldVel() * (p1E - PDat.getDeltaKE());

  sysMom += PDat.getDeltaP();
  
  if (PDat.getSpecies().getID() == species1)
    constDelGsp1 += PDat.getDeltaP();
}

void 
OPThermalDiffusionE::eventUpdate(const GlobalEvent& iEvent, 
				  const NEventData& PDat) 
{
  stream(iEvent.getdt());
  delG += impulseDelG(PDat);
  updateConstDelG(PDat);
}

void 
OPThermalDiffusionE::eventUpdate(const LocalEvent& iEvent, 
				  const NEventData& PDat) 
{
  stream(iEvent.getdt());
  delG += impulseDelG(PDat);
  updateConstDelG(PDat);
}

void 
OPThermalDiffusionE::eventUpdate(const System&, 
				  const NEventData& PDat, 
				  const double& edt) 
{ 
  stream(edt);
  delG += impulseDelG(PDat);
  updateConstDelG(PDat);
}
  
void 
OPThermalDiffusionE::eventUpdate(const IntEvent& iEvent, 
				  const PairEventData& PDat)
{
  stream(iEvent.getdt());
  delG += impulseDelG(PDat);
  updateConstDelG(PDat);
}

Vector  
OPThermalDiffusionE::impulseDelG(const NEventData& ndat) 
{ 
  Vector  acc(0,0,0);
    
  BOOST_FOREACH(const PairEventData& dat, ndat.L2partChanges)
    acc += impulseDelG(dat);
  
  return acc;
}

void 
OPThermalDiffusionE::updateConstDelG(const NEventData& ndat)
{
  BOOST_FOREACH(const ParticleEventData& dat, ndat.L1partChanges)
    updateConstDelG(dat);
  
  BOOST_FOREACH(const PairEventData& dat, ndat.L2partChanges)
    updateConstDelG(dat);
}
