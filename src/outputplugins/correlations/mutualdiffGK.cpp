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

#include "mutualdiffGK.hpp"
#include "../../dynamics/include.hpp"
#include "../1partproperty/kenergy.hpp"
#include "../0partproperty/misc.hpp"
#include "../../datatypes/vector.xml.hpp"

OPMutualDiffusionGK::OPMutualDiffusionGK(const DYNAMO::SimData* tmp, const XMLNode& XML):
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
OPMutualDiffusionGK::operator<<(const XMLNode& XML)
{
  try 
    {
      if (XML.isAttributeSet("Length"))
	CorrelatorLength = boost::lexical_cast<unsigned int>(XML.getAttribute("Length"));

      if (XML.isAttributeSet("dt"))
	dt = Sim->dynamics.units().unitTime() * 
	  boost::lexical_cast<Iflt>(XML.getAttribute("dt"));

	if (XML.isAttributeSet("t"))
	  dt = Sim->dynamics.units().unitTime() * 
	    boost::lexical_cast<Iflt>(XML.getAttribute("t")) / CorrelatorLength;
    }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in OPMutualDiffusionGK";
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
	  D_throw() << "You must set Species1 and Species2 for mutal diffusion\n"
		    << nex.what();
	}
    }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in OPMutualDiffusionGK";
    }
}

void 
OPMutualDiffusionGK::stream(const Iflt edt)
{
  //Now test if we've gone over the step time
  if (currentdt + edt >= dt)
    {
      newG();
      currentdt += edt - dt;
      
      while (currentdt >= dt)
	{
	  currentdt -= dt;
	  newG();
	}
    }
  else
    currentdt += edt;    
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
OPMutualDiffusionGK::eventUpdate(const System&, const NEventData& PDat, const Iflt& edt) 
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

Iflt 
OPMutualDiffusionGK::rescaleFactor()
{
  return 1.0 / (Sim->dynamics.units().unitMutualDiffusion()
		* count * Sim->dynamics.units().simVolume()
		* Sim->getOutputPlugin<OPKEnergy>()->getAvgkT());
}

void 
OPMutualDiffusionGK::output(xmlw::XmlStream& XML)
{
  Iflt factor = rescaleFactor();
  
  Vector  acc = 0.5*(accG[0] + accG[accG.size()-1]);

  for (unsigned int i = 1; i < accG.size()-1; i++)
    acc += accG[i];

  acc *= factor * dt / Sim->dynamics.units().unitTime();

  XML << xmlw::tag("Correlator")
      << xmlw::attr("name") << name
      << xmlw::attr("size") << accG.size()
      << xmlw::attr("dt") << dt / Sim->dynamics.units().unitTime()
      << xmlw::attr("LengthInMFT") << dt * accG.size() 
    / Sim->getOutputPlugin<OPMisc>()->getMFT()
      << xmlw::attr("simFactor") << factor
      << xmlw::attr("SampleCount") << count
      << xmlw::tag("Integral") << acc
      << xmlw::endtag("Integral")
      << xmlw::chardata();
    
  //GK correlators start at 0
  for (size_t i = 0; i < accG.size(); i++)
    {
      XML << i * dt / Sim->dynamics.units().unitTime();
      for (size_t j = 0; j < NDIM; j++)
	XML << "\t" << accG[i][j] * factor;
      XML << "\n";
    }
  
  
  XML << xmlw::endtag("Correlator"); 
}

void 
OPMutualDiffusionGK::initialise()
{
  species1 = Sim->dynamics.getSpecies(species1name).getID();
  species2 = Sim->dynamics.getSpecies(species2name).getID();

  Sim->getOutputPlugin<OPKEnergy>();
  Sim->getOutputPlugin<OPMisc>();
  
  accG.resize(CorrelatorLength, Vector  (0,0,0));
  G.resize(CorrelatorLength, Vector (0,0,0));
  dt = getdt();
  
  Iflt sysMass = 0.0;

  BOOST_FOREACH(const ClonePtr<Species>& sp, Sim->dynamics.getSpecies())
    sysMass += sp->getMass() * sp->getCount();
  
  BOOST_FOREACH(const Particle& part, Sim->particleList)
    {
      sysMom += part.getVelocity() * Sim->dynamics.getSpecies(part).getMass();
      
      if (Sim->dynamics.getSpecies()[species1]->isSpecies(part))
	delGsp1 += part.getVelocity();
      
      if (Sim->dynamics.getSpecies()[species2]->isSpecies(part))
	delGsp2 += part.getVelocity();
    }
  
  delGsp1 *= Sim->dynamics.getSpecies()[species1]->getMass();
  delGsp2 *= Sim->dynamics.getSpecies()[species2]->getMass();
  
  massFracSp1 = (Sim->dynamics.getSpecies()[species1]->getCount() 
		 * Sim->dynamics.getSpecies()[species1]->getMass()) 
    / sysMass; 

  massFracSp2 = (Sim->dynamics.getSpecies()[species2]->getCount() 
		 * Sim->dynamics.getSpecies()[species2]->getMass()) / sysMass;

  I_cout() << "dt set to " << dt / Sim->dynamics.units().unitTime();
}

std::list<Vector  > 
OPMutualDiffusionGK::getAvgAcc() const
{
  std::list<Vector  > tmp;
  
  BOOST_FOREACH(const Vector  &val, accG)
    tmp.push_back(val/((Iflt) count));
  
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

  accPass ();
}

void 
OPMutualDiffusionGK::accPass()
{
  ++count;
  
  for (size_t i = 0; i < CorrelatorLength; ++i)
    for (size_t j = 0; j < NDIM; ++j)
      accG[i][j] += (delGsp1[j] - (massFracSp1 * sysMom[j])) * (G[i][j] - (massFracSp2 * sysMom[j]));
}

Iflt 
OPMutualDiffusionGK::getdt()
{
  //Get the simulation temperature
    if (dt == 0.0)
      {
	if (Sim->lastRunMFT != 0.0)
	  return Sim->lastRunMFT * 50.0 / CorrelatorLength;
	else
	  return 5.0 / (((Iflt) CorrelatorLength)*sqrt(Sim->dynamics.getLiouvillean().getkT()) * CorrelatorLength);
      }
    else 
      return dt;
}

