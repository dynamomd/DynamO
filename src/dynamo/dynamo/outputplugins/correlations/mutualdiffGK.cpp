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

#include <dynamo/outputplugins/correlations/mutualdiffGK.hpp>
#include <dynamo/include.hpp>
#include <dynamo/outputplugins/1partproperty/kenergy.hpp>
#include <dynamo/outputplugins/0partproperty/misc.hpp>

namespace dynamo {
  OPMutualDiffusionGK::OPMutualDiffusionGK(const dynamo::Simulation* tmp, 
					   const magnet::xml::Node& XML):
    OutputPlugin(tmp, "MutualDiffusionGK", 60), //Note the sort order set later
    count(0),
    dt(0),
    currentdt(0.0),
    delGsp1(0,0,0),
    delGsp2(0,0,0),
    sysMom(0,0,0),
    massFracSp1(1),
    massFracSp2(1),
    CorrelatorLength(100),
    currCorrLen(0),
    notReady(true)

  {
    operator<<(XML);
  }

  void 
  OPMutualDiffusionGK::operator<<(const magnet::xml::Node& XML)
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
	M_throw() << "Failed a lexical cast in OPMutualDiffusionGK";
      }
  
    try 
      {
	try {
	  species1name = boost::lexical_cast<std::string>
	    (XML.getAttribute("Species1"));

	  species2name = boost::lexical_cast<std::string>
	    (XML.getAttribute("Species2"));

	} catch (std::exception& nex)
	  {
	    M_throw() << "You must set Species1 and Species2 for mutal diffusion\n"
		      << nex.what();
	  }
      }
    catch (boost::bad_lexical_cast &)
      {
	M_throw() << "Failed a lexical cast in OPMutualDiffusionGK";
      }
  }

  void 
  OPMutualDiffusionGK::stream(const double edt)
  {
    //Now test if we've gone over the step time
    currentdt += edt;
    while (currentdt >= dt)
      {
	newG();
	currentdt -= dt;
      }
  }

  void 
  OPMutualDiffusionGK::eventUpdate(const GlobalEvent& iEvent, const NEventData& PDat) 
  {
    stream(iEvent.getdt());
    updateDelG(PDat);
  }

  void 
  OPMutualDiffusionGK::eventUpdate(const LocalEvent& iEvent, const NEventData& PDat) 
  {
    stream(iEvent.getdt());
    updateDelG(PDat);
  }

  void 
  OPMutualDiffusionGK::eventUpdate(const System&, const NEventData& PDat, const double& edt) 
  { 
    stream(edt);
    updateDelG(PDat);
  }

  void 
  OPMutualDiffusionGK::eventUpdate(const IntEvent& iEvent, const PairEventData& PDat)
  {
    stream(iEvent.getdt());
    updateDelG(PDat);
  }

  double 
  OPMutualDiffusionGK::rescaleFactor()
  {
    return 1.0 / (Sim->units.unitMutualDiffusion()
		  * count * Sim->getSimVolume()
		  * Sim->getOutputPlugin<OPKEnergy>()->getAvgkT());
  }

  void 
  OPMutualDiffusionGK::output(magnet::xml::XmlStream& XML)
  {
    double factor = rescaleFactor();
  
    Vector  acc = 0.5*(accG[0] + accG[accG.size()-1]);

    for (unsigned int i = 1; i < accG.size()-1; i++)
      acc += accG[i];

    acc *= factor * dt / Sim->units.unitTime();

    XML << magnet::xml::tag("Correlator")
	<< magnet::xml::attr("name") << name
	<< magnet::xml::attr("size") << accG.size()
	<< magnet::xml::attr("dt") << dt / Sim->units.unitTime()
	<< magnet::xml::attr("LengthInMFT") << dt * accG.size() 
      / Sim->getOutputPlugin<OPMisc>()->getMFT()
	<< magnet::xml::attr("simFactor") << factor
	<< magnet::xml::attr("SampleCount") << count
	<< magnet::xml::tag("Integral") << acc
	<< magnet::xml::endtag("Integral")
	<< magnet::xml::chardata();
    
    //GK correlators start at 0
    for (size_t i = 0; i < accG.size(); i++)
      {
	XML << i * dt / Sim->units.unitTime();
	for (size_t j = 0; j < NDIM; j++)
	  XML << "\t" << accG[i][j] * factor;
	XML << "\n";
      }
  
  
    XML << magnet::xml::endtag("Correlator"); 
  }

  void 
  OPMutualDiffusionGK::initialise()
  {
    species1 = Sim->species[species1name].getID();
    species2 = Sim->species[species2name].getID();

    if (!(Sim->getOutputPlugin<OPMisc>()))
      M_throw() << "MutualDiffusionGK requires Misc output plugin!";
    if (!(Sim->getOutputPlugin<OPKEnergy>()))
      M_throw() << "MutualDiffusionGK requires KEnergy output plugin!";
  
    accG.resize(CorrelatorLength, Vector  (0,0,0));
    G.resize(CorrelatorLength, Vector (0,0,0));
    dt = getdt();
  
    double sysMass = 0.0;
  
    massFracSp1 = 0;
    massFracSp2 = 0;

    BOOST_FOREACH(const Particle& part, Sim->particleList)
      {
	double mass = Sim->species[part].getMass(part.getID());
	sysMom += part.getVelocity() * mass;
	sysMass += mass;

	if (Sim->species[species1]->isSpecies(part))
	  {
	    delGsp1 += part.getVelocity() * mass;
	    massFracSp1 += mass;
	  }
      
	if (Sim->species[species2]->isSpecies(part))
	  {
	    delGsp2 += part.getVelocity() * mass;
	    massFracSp2 += mass;
	  }
      }
  
    massFracSp1 /= sysMass; 
    massFracSp2 /= sysMass; 

    dout << "dt set to " << dt / Sim->units.unitTime() << std::endl;
  }

  std::list<Vector  > 
  OPMutualDiffusionGK::getAvgAcc() const
  {
    std::list<Vector  > tmp;
  
    BOOST_FOREACH(const Vector  &val, accG)
      tmp.push_back(val/((double) count));
  
    return tmp;
  }

  void 
  OPMutualDiffusionGK::updateDelG(const PairEventData& PDat) 
  {
    updateDelG(PDat.particle1_);
    updateDelG(PDat.particle2_);
  }

  void 
  OPMutualDiffusionGK::updateDelG(const ParticleEventData& PDat) 
  {
    sysMom += PDat.getDeltaP();
  
    if (PDat.getSpecies().getID() == species1)
      delGsp1 += PDat.getDeltaP();
  
    if (PDat.getSpecies().getID() == species2)
      delGsp2 += PDat.getDeltaP();
  
  }

  void 
  OPMutualDiffusionGK::updateDelG(const NEventData& ndat)
  {
    BOOST_FOREACH(const ParticleEventData& dat, ndat.L1partChanges)
      updateDelG(dat);
  
    BOOST_FOREACH(const PairEventData& dat, ndat.L2partChanges)
      updateDelG(dat);
  }

  void 
  OPMutualDiffusionGK::newG()
  {
    G.push_front(delGsp2);

    //This ensures the list gets to accumilator size
    if (notReady)
      {
	if (++currCorrLen != CorrelatorLength)
	  return;
      
	notReady = false;
      }

    accPass();
  }

  void 
  OPMutualDiffusionGK::accPass()
  {
    ++count;
  
    for (size_t i = 0; i < CorrelatorLength; ++i)
      for (size_t j = 0; j < NDIM; ++j)
	accG[i][j] += (delGsp1[j] - (massFracSp1 * sysMom[j])) * (G[i][j] - (massFracSp2 * sysMom[j]));
  }

  double 
  OPMutualDiffusionGK::getdt()
  {
    //Get the simulation temperature
    if (dt == 0.0)
      {
	if (Sim->lastRunMFT != 0.0)
	  return Sim->lastRunMFT * 50.0 / CorrelatorLength;
	else
	  return 5.0 / (((double) CorrelatorLength)*sqrt(Sim->dynamics->getkT()) * CorrelatorLength);
      }
    else 
      return dt;
  }
}
