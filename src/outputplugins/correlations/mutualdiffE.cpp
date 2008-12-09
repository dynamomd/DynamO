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

#include "mutualdiffE.hpp"
#include "../../dynamics/include.hpp"
#include "../1partproperty/kenergy.hpp"
#include "../0partproperty/misc.hpp"
#include "../../datatypes/vector.xml.hpp"

COPMutualDiffusionE::COPMutualDiffusionE(const DYNAMO::SimData* tmp, const XMLNode& XML):
  COutputPlugin(tmp, "MutualDiffusionE", 60), //Note the sort order set later
  count(0),
  dt(0),
  currentdt(0.0),
  delGsp1(0.0),
  delGsp2(0.0),
  Gsp1(0.0),
  Gsp2(0.0),
  sysMom(0.0),
  massFracSp1(1),
  massFracSp2(1),
  CorrelatorLength(100),
  currCorrLen(0),
  notReady(true)

{
  operator<<(XML);
}

void 
COPMutualDiffusionE::operator<<(const XMLNode& XML)
{
  try 
    {
      if (XML.isAttributeSet("Length"))
	CorrelatorLength = boost::lexical_cast<unsigned int>(XML.getAttribute("Length"));

      if (XML.isAttributeSet("dt"))
	dt = Sim->Dynamics.units().unitTime() * 
	  boost::lexical_cast<Iflt>(XML.getAttribute("dt"));

	if (XML.isAttributeSet("t"))
	  dt = Sim->Dynamics.units().unitTime() * 
	    boost::lexical_cast<Iflt>(XML.getAttribute("t")) / CorrelatorLength;
    }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in COPMutualDiffusionE";
    }
  
  try 
    {
      try {
	species1 = Sim->Dynamics
	  .getSpecies(boost::lexical_cast<std::string>
		      (XML.getAttribute("Species1"))).getID();

	species2 = Sim->Dynamics
	  .getSpecies(boost::lexical_cast<std::string>
		      (XML.getAttribute("Species2"))).getID();

      }
      catch (boost::bad_lexical_cast &)
	{
	  D_throw() << "Failed a lexical cast in COPMutualDiffusionE";
	}

    } catch (std::exception& nex)
    {
      D_throw() << "You must set Species1 and Species2 for mutal diffusion\n"
		<< nex.what();
    }
}

void 
COPMutualDiffusionE::stream(const Iflt edt)
{
  CVector<> grad1 = (delGsp1 - (massFracSp1 * sysMom)),
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
COPMutualDiffusionE::eventUpdate(const CGlobEvent& iEvent,
				 const CNParticleData& PDat) 
{
  stream(iEvent.getdt());
  updateDelG(PDat);
}

void 
COPMutualDiffusionE::eventUpdate(const CLocalEvent& iEvent, 
				 const CNParticleData& PDat) 
{
  stream(iEvent.getdt());
  updateDelG(PDat);
}

void 
COPMutualDiffusionE::eventUpdate(const CSystem&, const CNParticleData& PDat,
				 const Iflt& edt) 
{ 
  stream(edt);
  updateDelG(PDat);
}

void 
COPMutualDiffusionE::eventUpdate(const CIntEvent& iEvent, 
				 const C2ParticleData& PDat)
{
  stream(iEvent.getdt());
  updateDelG(PDat);
}

Iflt 
COPMutualDiffusionE::rescaleFactor()
{
  return 0.5 / (Sim->Dynamics.units().unitTime()
		* Sim->Dynamics.units().unitMutualDiffusion()
		* count * Sim->Dynamics.units().simVolume()
		* Sim->getOutputPlugin<COPKEnergy>()->getAvgkT());
}

void 
COPMutualDiffusionE::output(xmlw::XmlStream& XML)
{
  Iflt factor = rescaleFactor();

  XML << xmlw::tag("EinsteinCorrelator")
      << xmlw::attr("name") << name
      << xmlw::attr("size") << accG.size()
      << xmlw::attr("dt") << dt / Sim->Dynamics.units().unitTime()
      << xmlw::attr("LengthInMFT") << dt * accG.size() 
    / Sim->getOutputPlugin<COPMisc>()->getMFT()
      << xmlw::attr("simFactor") << factor
      << xmlw::attr("SampleCount") << count
      << xmlw::chardata();
    
  for (size_t i = 0; i < accG.size(); i++)
    {
      XML << (i + 1) * dt / Sim->Dynamics.units().unitTime();
      for (size_t j = 0; j < NDIM; j++)
	XML << "\t" << accG[i][j] * factor;
      XML << "\n";
    }
  
  
  XML << xmlw::endtag("EinsteinCorrelator"); 
}

void 
COPMutualDiffusionE::initialise()
{
  Sim->getOutputPlugin<COPKEnergy>();
  Sim->getOutputPlugin<COPMisc>();
  
  accG.resize(CorrelatorLength, CVector<> (0.0));

  G1.resize(CorrelatorLength, CVector<>(0.0));

  G2.resize(CorrelatorLength, CVector<>(0.0));

  dt = getdt();
  
  Iflt sysMass = 0.0;

  BOOST_FOREACH(const CSpecies& sp, Sim->Dynamics.getSpecies())
    sysMass += sp.getMass() * sp.getCount();
  
  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    {
      sysMom += part.getVelocity() * Sim->Dynamics.getSpecies(part).getMass();
      
      if (Sim->Dynamics.getSpecies()[species1].isSpecies(part))
	delGsp1 += part.getVelocity();
      
      if (Sim->Dynamics.getSpecies()[species2].isSpecies(part))
	delGsp2 += part.getVelocity();
    }
  
  delGsp1 *= Sim->Dynamics.getSpecies()[species1].getMass();
  delGsp2 *= Sim->Dynamics.getSpecies()[species2].getMass();
  
  massFracSp1 = (Sim->Dynamics.getSpecies()[species1].getCount() 
		 * Sim->Dynamics.getSpecies()[species1].getMass()) 
    / sysMass; 

  massFracSp2 = (Sim->Dynamics.getSpecies()[species2].getCount() 
		 * Sim->Dynamics.getSpecies()[species2].getMass()) / sysMass;

  I_cout() << "dt set to " << dt / Sim->Dynamics.units().unitTime();
}

std::list<CVector<> > 
COPMutualDiffusionE::getAvgAcc() const
{
  std::list<CVector<> > tmp;
  
  BOOST_FOREACH(const CVector<> &val, accG)
    tmp.push_back(val/((Iflt) count));
  
  return tmp;
}

void 
COPMutualDiffusionE::updateDelG(const C2ParticleData& PDat) 
{
  updateDelG(PDat.particle1_);
  updateDelG(PDat.particle2_);
}

void 
COPMutualDiffusionE::updateDelG(const C1ParticleData& PDat) 
{
  sysMom += PDat.getDeltaP();
  
  if (PDat.getSpecies().getID() == species1)
    delGsp1 += PDat.getDeltaP();
  
  if (PDat.getSpecies().getID() == species2)
    delGsp2 += PDat.getDeltaP();
  
}

void 
COPMutualDiffusionE::updateDelG(const CNParticleData& ndat)
{
  BOOST_FOREACH(const C1ParticleData& dat, ndat.L1partChanges)
    updateDelG(dat);
  
  BOOST_FOREACH(const C2ParticleData& dat, ndat.L2partChanges)
    updateDelG(dat);
}

void 
COPMutualDiffusionE::newG()
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
COPMutualDiffusionE::accPass()
{
  ++count;
  
  CVector<> sumSp1(0), sumSp2(0);

  for (size_t i = 0; i < CorrelatorLength; ++i)
    {
      sumSp1 += G1[i];
      sumSp2 += G2[i];

      accG[i] += sumSp1 * sumSp2;
    }
}

Iflt 
COPMutualDiffusionE::getdt()
{
  //Get the simulation temperature
    if (dt == 0.0)
      {
	if (Sim->lastRunMFT != 0.0)
	  return Sim->lastRunMFT * 50.0 / CorrelatorLength;
	else
	  return 5.0 / (((Iflt) CorrelatorLength)*sqrt(Sim->Dynamics.getkT()) * CorrelatorLength);
      }
    else 
      return dt;
}

