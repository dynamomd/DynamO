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

#include <dynamo/outputplugins/correlations/thermalCondSpeciesSpeciesE.hpp>
#include <dynamo/dynamics/include.hpp>
#include <dynamo/dynamics/interactions/intEvent.hpp>
#include <dynamo/outputplugins/1partproperty/kenergy.hpp>
#include <dynamo/simulation/ensemble.hpp>
#include <dynamo/outputplugins/0partproperty/misc.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>

namespace dynamo {
  OPThermalConductivitySpeciesSpeciesE::OPThermalConductivitySpeciesSpeciesE(const dynamo::SimData* tmp,
									     const magnet::xml::Node& XML):
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
  OPThermalConductivitySpeciesSpeciesE::operator<<(const magnet::xml::Node& XML)
  {
    try 
      {
	if (XML.hasAttribute("Length"))
	  CorrelatorLength = XML.getAttribute("Length").as<size_t>();
      
	if (XML.hasAttribute("dt"))
	  dt = Sim->dynamics.units().unitTime() * 
	    XML.getAttribute("dt").as<double>();
      
	if (XML.hasAttribute("t"))
	  dt = Sim->dynamics.units().unitTime() * 
	    XML.getAttribute("t").as<double>() / CorrelatorLength;
      }
    catch (boost::bad_lexical_cast &)
      {
	M_throw() << "Failed a lexical cast in OPCorrelator";
      }
  
  }

  void 
  OPThermalConductivitySpeciesSpeciesE::initialise()
  {
    size_t Nsp(Sim->species.size());
  
    constDelG.resize(Nsp, Vector (0,0,0));
  
    delG.resize(Nsp, Vector (0,0,0));
  
    G.resize(Nsp, boost::circular_buffer<Vector  >(CorrelatorLength, Vector(0,0,0)));
  
    accG2.resize(Nsp * Nsp, std::vector<Vector  >(CorrelatorLength, Vector(0,0,0)));
  
    if (!(Sim->getOutputPlugin<OPMisc>()))
      M_throw() << "ThermalConductivitySpeciesSpeciesE requires Misc output plugin!";
    if (!(Sim->getOutputPlugin<OPKEnergy>()))
      M_throw() << "ThermalConductivitySpeciesSpeciesE requires KEnergy output plugin!";
  
    if (dynamic_cast<const dynamo::EnsembleNVE* >(Sim->ensemble.get()) == NULL)
      M_throw() << "WARNING: This is only valid in the microcanonical"
	" ensemble!\nSee J.J. Erpenbeck, Phys. Rev. A 39, 4718 (1989) for more"
	"\n Essentially you need entropic data too for other ensembles";

    if (dt == 0.0)
      {
	if (Sim->lastRunMFT != 0.0)
	  dt = Sim->lastRunMFT * 50.0 / CorrelatorLength;
	else
	  dt = 10.0 / (((double) CorrelatorLength) 
		       * sqrt(Sim->liouvillean->getkT()) * CorrelatorLength);
      }
  
    //Sum up the constant Del G.
    BOOST_FOREACH(const shared_ptr<Species>& spec, Sim->species)
      BOOST_FOREACH(const size_t& id, *spec->getRange())
      {
	const Particle& part(Sim->particleList[id]);
	constDelG[spec->getID()] += part.getVelocity() 
	  * Sim->liouvillean->getParticleKineticEnergy(part);
      }
  
    dout << "dt set to " << dt / Sim->dynamics.units().unitTime() << std::endl;
  }

  double 
  OPThermalConductivitySpeciesSpeciesE::rescaleFactor() 
  { 
    return Sim->dynamics.units().unitk() 
      /(//This next line should be 1 however we have scaled the
	//correlator time as well
	Sim->dynamics.units().unitTime() 
	* Sim->dynamics.units().unitThermalCond() * 2.0 
	* count * pow(Sim->getOutputPlugin<OPKEnergy>()->getAvgkT(), 2)
	* Sim->dynamics.getSimVolume());
  }

  void 
  OPThermalConductivitySpeciesSpeciesE::output(magnet::xml::XmlStream &XML)
  {
    XML << magnet::xml::tag("EinsteinCorrelator")
	<< magnet::xml::attr("name") << name
	<< magnet::xml::attr("size") << accG2.size()
	<< magnet::xml::attr("dt") << dt/Sim->dynamics.units().unitTime()
	<< magnet::xml::attr("LengthInMFT") << dt * accG2.size()
      / Sim->getOutputPlugin<OPMisc>()->getMFT()
	<< magnet::xml::attr("simFactor") << rescaleFactor()
	<< magnet::xml::attr("SampleCount") << count;
  
    double factor = rescaleFactor();

    size_t Nsp(Sim->species.size());
  
    for (size_t id1(0); id1 < Nsp; ++id1)
      for (size_t id2(0); id2 < Nsp; ++id2)      
	{
	  XML << magnet::xml::tag("Component") 
	      << magnet::xml::attr("Species1") << id1
	      << magnet::xml::attr("Species2") << id2
	      << magnet::xml::chardata();

	  for (unsigned int i = 0; i < accG2.size(); i++)
	    {
	      XML << (1+i) * dt / Sim->dynamics.units().unitTime()
		  << "\t ";
	    
	      for (size_t j=0;j<NDIM;j++)
		XML << accG2[id1 + Nsp * id2][i][j] * factor 
		    << "\t ";
	    
	      XML << "\n";
	    }

	  XML << magnet::xml::endtag("Component");
	}
  
    XML << magnet::xml::endtag("EinsteinCorrelator");
  }

  void 
  OPThermalConductivitySpeciesSpeciesE::stream(const double& edt)
  {
    size_t Nsp(Sim->species.size());
  
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
  OPThermalConductivitySpeciesSpeciesE::eventUpdate(const GlobalEvent& iEvent, 
						    const NEventData& PDat) 
  {
    stream(iEvent.getdt());
    //impulseDelG(PDat);
    updateConstDelG(PDat);
  }

  void 
  OPThermalConductivitySpeciesSpeciesE::eventUpdate(const LocalEvent& iEvent, 
						    const NEventData& PDat) 
  {
    stream(iEvent.getdt());
    //impulseDelG(PDat);
    updateConstDelG(PDat);
  }

  void 
  OPThermalConductivitySpeciesSpeciesE::eventUpdate(const System&, 
						    const NEventData& PDat,
						    const double& edt) 
  { 
    stream(edt);
    //impulseDelG(PDat);
    updateConstDelG(PDat);
  }
  
  void 
  OPThermalConductivitySpeciesSpeciesE::eventUpdate(const IntEvent& iEvent,
						    const PairEventData& PDat)
  {
    stream(iEvent.getdt());
    //impulseDelG(PDat);
    updateConstDelG(PDat);
  }

  void
  OPThermalConductivitySpeciesSpeciesE::impulseDelG(const NEventData& ndat) 
  { 
    /*  Vector  acc(0);
  
	BOOST_FOREACH(const C2ParticleData& dat, ndat.L2partChanges)
	acc += impulseDelG(dat);
  
	return acc;*/
  }

  void
  OPThermalConductivitySpeciesSpeciesE::impulseDelG(const PairEventData& PDat)
  {
    //return PDat.rij * PDat.particle1_.getDeltaeCalc();
  }

  void 
  OPThermalConductivitySpeciesSpeciesE::newG()
  {
    //This ensures the list stays at accumilator size
  
    size_t Nsp(Sim->species.size());
  
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

    size_t Nsp(Sim->species.size());
  
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
  OPThermalConductivitySpeciesSpeciesE::updateConstDelG(const NEventData& ndat)
  {
    BOOST_FOREACH(const ParticleEventData& dat, ndat.L1partChanges)
      updateConstDelG(dat);
  
    BOOST_FOREACH(const PairEventData& dat, ndat.L2partChanges)
      updateConstDelG(dat);
  }

  void 
  OPThermalConductivitySpeciesSpeciesE::updateConstDelG(const ParticleEventData& PDat)
  {
    double p1E = Sim->liouvillean->getParticleKineticEnergy(PDat.getParticle());
  
    constDelG[PDat.getSpecies().getID()] += PDat.getParticle().getVelocity() * p1E 
      - PDat.getOldVel() * (p1E - PDat.getDeltaKE());
  }

  void 
  OPThermalConductivitySpeciesSpeciesE::updateConstDelG(const PairEventData& PDat)
  {
    double p1E = Sim->liouvillean->getParticleKineticEnergy(PDat.particle1_.getParticle());
    double p2E = Sim->liouvillean->getParticleKineticEnergy(PDat.particle2_.getParticle());
  
    constDelG[Sim->species[PDat.particle1_.getParticle()].getID()]
      += PDat.particle1_.getParticle().getVelocity() * p1E
      - PDat.particle1_.getOldVel() * (p1E - PDat.particle1_.getDeltaKE());
  
    constDelG[Sim->species[PDat.particle2_.getParticle()].getID()]
      += PDat.particle2_.getParticle().getVelocity() * p2E
      - PDat.particle2_.getOldVel() * (p2E - PDat.particle2_.getDeltaKE());
  }
}
