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

#include <dynamo/outputplugins/correlations/mutualdiffE.hpp>
#include <dynamo/dynamics/include.hpp>
#include <dynamo/outputplugins/1partproperty/kenergy.hpp>
#include <dynamo/outputplugins/0partproperty/misc.hpp>

namespace dynamo {
  OPMutualDiffusionE::OPMutualDiffusionE(const dynamo::SimData* tmp, 
					 const magnet::xml::Node& XML):
    OutputPlugin(tmp, "MutualDiffusionE", 60), //Note the sort order set later
    count(0),
    dt(0.0),
    currentdt(0.0),
    delGsp1(0,0,0),
    delGsp2(0,0,0),
    Gsp1(0,0,0),
    Gsp2(0,0,0),
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
  OPMutualDiffusionE::operator<<(const magnet::xml::Node& XML)
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
	M_throw() << "Failed a lexical cast in OPMutualDiffusionE";
      }
  
    try 
      {
	try {
	  species1name = boost::lexical_cast<std::string>
	    (XML.getAttribute("Species1"));

	  species2name = boost::lexical_cast<std::string>
	    (XML.getAttribute("Species2"));

	}
	catch (boost::bad_lexical_cast &)
	  {
	    M_throw() << "Failed a lexical cast in OPMutualDiffusionE";
	  }

      } catch (std::exception& nex)
      {
	M_throw() << "You must set Species1 and Species2 for mutal diffusion\n"
		  << nex.what();
      }
  }

  void 
  OPMutualDiffusionE::stream(const double edt)
  {
    Vector  grad1 = (delGsp1 - (massFracSp1 * sysMom)),
      grad2 = (delGsp2 - (massFracSp2 * sysMom));  

    //Now test if we've gone over the step time
    if (currentdt + edt >= dt)
      {
	Gsp1 += grad1 * (dt - currentdt);
	Gsp2 += grad2 * (dt - currentdt);

	newG();

	currentdt += edt - dt;
      
	while (currentdt >= dt)
	  {
	    Gsp1 = grad1 * dt;	  
	    Gsp2 = grad2 * dt;

	    currentdt -= dt;

	    newG();
	  }

	Gsp1 = grad1 * currentdt;    
	Gsp2 = grad2 * currentdt;
      }
    else
      {
	Gsp1 += grad1 * edt;      
	Gsp2 += grad2 * edt;

	currentdt += edt;    
      }
  }

  void 
  OPMutualDiffusionE::eventUpdate(const GlobalEvent& iEvent,
				  const NEventData& PDat) 
  {
    stream(iEvent.getdt());
    updateDelG(PDat);
  }

  void 
  OPMutualDiffusionE::eventUpdate(const LocalEvent& iEvent, 
				  const NEventData& PDat) 
  {
    stream(iEvent.getdt());
    updateDelG(PDat);
  }

  void 
  OPMutualDiffusionE::eventUpdate(const System&, const NEventData& PDat,
				  const double& edt) 
  { 
    stream(edt);
    updateDelG(PDat);
  }

  void 
  OPMutualDiffusionE::eventUpdate(const IntEvent& iEvent, 
				  const PairEventData& PDat)
  {
    stream(iEvent.getdt());
    updateDelG(PDat);
  }

  double 
  OPMutualDiffusionE::rescaleFactor()
  {
    return 0.5 / (Sim->dynamics.units().unitTime()
		  * Sim->dynamics.units().unitMutualDiffusion()
		  * count * Sim->dynamics.getSimVolume()
		  * Sim->getOutputPlugin<OPKEnergy>()->getAvgkT());
  }

  void 
  OPMutualDiffusionE::output(magnet::xml::XmlStream& XML)
  {
    double factor = rescaleFactor();

    XML << magnet::xml::tag("EinsteinCorrelator")
	<< magnet::xml::attr("name") << name
	<< magnet::xml::attr("size") << accG.size()
	<< magnet::xml::attr("dt") << dt / Sim->dynamics.units().unitTime()
	<< magnet::xml::attr("LengthInMFT") << dt * accG.size() 
      / Sim->getOutputPlugin<OPMisc>()->getMFT()
	<< magnet::xml::attr("simFactor") << factor
	<< magnet::xml::attr("SampleCount") << count
	<< magnet::xml::chardata();
    
    for (size_t i = 0; i < accG.size(); i++)
      {
	XML << (i + 1) * dt / Sim->dynamics.units().unitTime();
	for (size_t j = 0; j < NDIM; j++)
	  XML << "\t" << accG[i][j] * factor;
	XML << "\n";
      }
  
  
    XML << magnet::xml::endtag("EinsteinCorrelator"); 
  }

  void 
  OPMutualDiffusionE::initialise()
  {
    species1 = Sim->dynamics.getSpecies(species1name).getID();
    species2 = Sim->dynamics.getSpecies(species2name).getID();
  
    Sim->getOutputPlugin<OPKEnergy>();
    Sim->getOutputPlugin<OPMisc>();
  
    accG.resize(CorrelatorLength, Vector  (0,0,0));

    G1.resize(CorrelatorLength, Vector (0,0,0));

    G2.resize(CorrelatorLength, Vector (0,0,0));

    dt = getdt();
  
    double sysMass = 0.0;
  
    massFracSp1 = 0;
    massFracSp2 = 0;

    BOOST_FOREACH(const Particle& part, Sim->particleList)
      {
	double mass = Sim->dynamics.getSpecies(part).getMass(part.getID());
	sysMass += mass;
	sysMom += part.getVelocity() * mass;
      
	if (Sim->dynamics.getSpecies()[species1]->isSpecies(part))
	  {
	    delGsp1 += part.getVelocity() * mass;
	    massFracSp1 += mass;
	  }
      
	if (Sim->dynamics.getSpecies()[species2]->isSpecies(part))
	  {
	    delGsp2 += part.getVelocity() * mass;
	    massFracSp2 += mass;
	  }
      }
    
    massFracSp1 /= sysMass; 
    massFracSp2 /= sysMass;

    dout << "dt set to " << dt / Sim->dynamics.units().unitTime() << std::endl;
  }

  std::list<Vector  > 
  OPMutualDiffusionE::getAvgAcc() const
  {
    std::list<Vector  > tmp;
  
    BOOST_FOREACH(const Vector  &val, accG)
      tmp.push_back(val/((double) count));
  
    return tmp;
  }

  void 
  OPMutualDiffusionE::updateDelG(const PairEventData& PDat) 
  {
    updateDelG(PDat.particle1_);
    updateDelG(PDat.particle2_);
  }

  void 
  OPMutualDiffusionE::updateDelG(const ParticleEventData& PDat) 
  {
    sysMom += PDat.getDeltaP();
  
    if (PDat.getSpecies().getID() == species1)
      delGsp1 += PDat.getDeltaP();
  
    if (PDat.getSpecies().getID() == species2)
      delGsp2 += PDat.getDeltaP();
  
  }

  void 
  OPMutualDiffusionE::updateDelG(const NEventData& ndat)
  {
    BOOST_FOREACH(const ParticleEventData& dat, ndat.L1partChanges)
      updateDelG(dat);
  
    BOOST_FOREACH(const PairEventData& dat, ndat.L2partChanges)
      updateDelG(dat);
  }

  void 
  OPMutualDiffusionE::newG()
  {
    G1.push_front(Gsp1);
    G2.push_front(Gsp2);

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
  OPMutualDiffusionE::accPass()
  {
    ++count;
  
    Vector  sumSp1(0,0,0), sumSp2(0,0,0);

    for (size_t i = 0; i < CorrelatorLength; ++i)
      {
	sumSp1 += G1[i];
	sumSp2 += G2[i];

	for (size_t j(0); j < NDIM; ++j)
	  accG[i][j] += sumSp1[j] * sumSp2[j];
      }
  }

  double 
  OPMutualDiffusionE::getdt()
  {
    //Get the simulation temperature
    if (dt == 0.0)
      {
	if (Sim->lastRunMFT != 0.0)
	  return Sim->lastRunMFT * 50.0 / CorrelatorLength;
	else
	  return 5.0 / (((double) CorrelatorLength)*sqrt(Sim->dynamics.getLiouvillean().getkT()) * CorrelatorLength);
      }
    else 
      return dt;
  }
}
