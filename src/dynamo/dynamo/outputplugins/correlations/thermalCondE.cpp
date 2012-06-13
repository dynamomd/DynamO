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

#include <dynamo/outputplugins/correlations/thermalCondE.hpp>
#include <dynamo/include.hpp>
#include <dynamo/interactions/intEvent.hpp>
#include <dynamo/outputplugins/1partproperty/kenergy.hpp>
#include <dynamo/ensemble.hpp>
#include <dynamo/outputplugins/0partproperty/misc.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>

namespace dynamo {
  OPThermalConductivityE::OPThermalConductivityE(const dynamo::SimData* tmp,
						 const magnet::xml::Node& XML):
    OutputPlugin(tmp,"ThermalConductivityE"),
    G(100),
    count(0),
    dt(0),
    currentdt(0.0),
    constDelG(0,0,0), 
    delG(0,0,0),
    currlen(0),
    notReady(true),
    CorrelatorLength(100)
  {
    operator<<(XML);
  }

  void 
  OPThermalConductivityE::operator<<(const magnet::xml::Node& XML)
  {
    try 
      {
	if (XML.hasAttribute("Length"))
	  CorrelatorLength = XML.getAttribute("Length").as<size_t>();
      
	if (XML.hasAttribute("dt"))
	  dt = Sim->units.unitTime() * 
	    XML.getAttribute("dt").as<double>();
      
	if (XML.hasAttribute("t"))
	  dt = Sim->units.unitTime() * 
	    XML.getAttribute("t").as<double>() / CorrelatorLength;
      }
    catch (boost::bad_lexical_cast &)
      {
	M_throw() << "Failed a lexical cast in OPCorrelator";
      }
  
  }

  void 
  OPThermalConductivityE::initialise()
  {
    G.resize(CorrelatorLength, Vector (0,0,0));
    accG2.resize(CorrelatorLength, Vector (0,0,0));
    if (!(Sim->getOutputPlugin<OPMisc>()))
      M_throw() << "ThermalConductivityE requires Misc output plugin!";
    if (!(Sim->getOutputPlugin<OPKEnergy>()))
      M_throw() << "ThermalConductivityE requires KEnergy output plugin!";
  
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
    BOOST_FOREACH(const Particle& part, Sim->particleList)
      constDelG += part.getVelocity () * Sim->liouvillean->getParticleKineticEnergy(part);
  
    dout << "dt set to " << dt / Sim->units.unitTime() << std::endl;
  }

  double 
  OPThermalConductivityE::rescaleFactor() 
  { 
    return Sim->units.unitk() 
      /(//This next line should be 1 however we have scaled the
	//correlator time as well
	Sim->units.unitTime() 
	* Sim->units.unitThermalCond() * 2.0 
	* count * pow(Sim->getOutputPlugin<OPKEnergy>()->getAvgkT(), 2)
	* Sim->getSimVolume());
  }

  void 
  OPThermalConductivityE::output(magnet::xml::XmlStream &XML)
  {
    XML << magnet::xml::tag("EinsteinCorrelator")
	<< magnet::xml::attr("name") << name
	<< magnet::xml::attr("size") << accG2.size()
	<< magnet::xml::attr("dt") << dt/Sim->units.unitTime()
	<< magnet::xml::attr("LengthInMFT") << dt * accG2.size()
      / Sim->getOutputPlugin<OPMisc>()->getMFT()
	<< magnet::xml::attr("simFactor") << rescaleFactor()
	<< magnet::xml::attr("SampleCount") << count
	<< magnet::xml::chardata();
  
    double factor = rescaleFactor();
  
    for (unsigned int i = 0; i < accG2.size(); i++)
      {
	XML   << (1+i) * dt / Sim->units.unitTime()
	      << "\t ";
      
	for (size_t j=0;j<NDIM;j++)
	  XML << accG2[i][j] * factor 
	      << "\t ";
      
	XML << "\n";
      }
  
    XML << magnet::xml::endtag("EinsteinCorrelator");
  }

  Vector  
  OPThermalConductivityE::impulseDelG(const PairEventData& PDat)
  {
    return PDat.rij * PDat.particle1_.getDeltaKE();
  }

  void 
  OPThermalConductivityE::updateConstDelG(const PairEventData& PDat)
  {
    double p1E = Sim->liouvillean->getParticleKineticEnergy(PDat.particle1_.getParticle());
    double p2E = Sim->liouvillean->getParticleKineticEnergy(PDat.particle2_.getParticle());
  
    constDelG += PDat.particle1_.getParticle().getVelocity() * p1E 
      + PDat.particle2_.getParticle().getVelocity() * p2E
      - PDat.particle1_.getOldVel() * (p1E - PDat.particle1_.getDeltaKE())
      - PDat.particle2_.getOldVel() * (p2E - PDat.particle2_.getDeltaKE());
  }

  void 
  OPThermalConductivityE::stream(const double& edt)
  {
    //Now test if we've gone over the step time
    if (currentdt + edt >= dt)
      {
	delG += constDelG * (dt - currentdt);
	newG();
	currentdt += edt - dt;

	while (currentdt >= dt)
	  {
	    delG = constDelG * dt;
	    currentdt -= dt;
	    newG();
	  }

	//Now calculate the start of the new delG
	delG = constDelG * currentdt;
      }
    else
      {
	currentdt += edt;
	delG += constDelG * edt;
      }
  }

  void 
  OPThermalConductivityE::eventUpdate(const GlobalEvent& iEvent, 
				      const NEventData& PDat) 
  {
    stream(iEvent.getdt());
    delG += impulseDelG(PDat);
    updateConstDelG(PDat);
  }

  void 
  OPThermalConductivityE::eventUpdate(const LocalEvent& iEvent, 
				      const NEventData& PDat) 
  {
    stream(iEvent.getdt());
    delG += impulseDelG(PDat);
    updateConstDelG(PDat);
  }

  void 
  OPThermalConductivityE::eventUpdate(const System&, 
				      const NEventData& PDat,
				      const double& edt) 
  { 
    stream(edt);
    delG += impulseDelG(PDat);
    updateConstDelG(PDat);
  }
  
  void 
  OPThermalConductivityE::eventUpdate(const IntEvent& iEvent,
				      const PairEventData& PDat)
  {
    stream(iEvent.getdt());
    delG += impulseDelG(PDat);
    updateConstDelG(PDat);
  }

  Vector  
  OPThermalConductivityE::impulseDelG(const NEventData& ndat) 
  { 
    Vector  acc(0,0,0);
  
    BOOST_FOREACH(const PairEventData& dat, ndat.L2partChanges)
      acc += impulseDelG(dat);
  
    return acc;
  }

  void 
  OPThermalConductivityE::newG()
  {
    //This ensures the list stays at accumilator size
  
    G.push_front(delG);
  
    if (notReady)
      {
	if (++currlen != CorrelatorLength)
	  return;
      
	notReady = false;
      }
  
    accPass();
  }

  void 
  OPThermalConductivityE::accPass()
  {
    ++count;
    Vector  sum(0,0,0);
  
    for (size_t i = 0; i < CorrelatorLength; ++i)
      {
	sum += G[i];

	for (size_t j(0); j < NDIM; ++j)
	  accG2[i][j] += sum[j] * sum[j];
      }
  }

  void 
  OPThermalConductivityE::updateConstDelG(const NEventData& ndat)
  {
    BOOST_FOREACH(const ParticleEventData& dat, ndat.L1partChanges)
      updateConstDelG(dat);
  
    BOOST_FOREACH(const PairEventData& dat, ndat.L2partChanges)
      updateConstDelG(dat);
  }

  void 
  OPThermalConductivityE::updateConstDelG(const ParticleEventData& PDat)
  {
    double p1E = Sim->liouvillean->getParticleKineticEnergy(PDat.getParticle());
  
    constDelG += PDat.getParticle().getVelocity() * p1E 
      - PDat.getOldVel() * (p1E - PDat.getDeltaKE());
  }
}
