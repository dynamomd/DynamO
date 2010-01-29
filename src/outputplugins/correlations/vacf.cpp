/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "vacf.hpp"
#include "../../dynamics/liouvillean/SLLOD.hpp"

OPVACF::OPVACF(const DYNAMO::SimData* tmp,const XMLNode& XML):
  OutputPlugin(tmp, "VACF", 60), //Note the sort order set later
  count(0),
  dt(0),
  currentdt(0.0),
  CorrelatorLength(100),
  currCorrLen(0),
  notReady(true)
{
  operator<<(XML);
}

void 
OPVACF::initialise()
{
  Sim->getOutputPlugin<OPMisc>();

  dt = getdt();

  G.resize(Sim->lN, boost::circular_buffer<Vector  >(CorrelatorLength, Vector(0,0,0)));

  accG2.resize(Sim->dynamics.getSpecies().size());

  BOOST_FOREACH(std::vector<Vector  >& listref, accG2)
    listref.resize(CorrelatorLength, Vector (0,0,0));

  I_cout() << "dt set to " << dt / Sim->dynamics.units().unitTime();
}

void 
OPVACF::operator<<(const XMLNode& XML)
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
      D_throw() << "Failed a lexical cast in OPVACF";
    }  
}

void 
OPVACF::eventUpdate(const CGlobEvent& iEvent, const CNParticleData& PDat) 
{
  //Move the time forward
  currentdt += iEvent.getdt();
  
  //Now test if we've gone over the step time
  while (currentdt >= dt)
    {
      currentdt -= dt;
      newG(PDat);
    }
}

void 
OPVACF::eventUpdate(const CLocalEvent& iEvent, const CNParticleData& PDat) 
{
  //Move the time forward
  currentdt += iEvent.getdt();
  
  //Now test if we've gone over the step time
  while (currentdt >= dt)
    {
      currentdt -= dt;
      newG(PDat);
    }
}

void 
OPVACF::eventUpdate(const CSystem&, const CNParticleData& PDat, const Iflt& edt) 
{ 
  //Move the time forward
  currentdt += edt;
  
  //Now test if we've gone over the step time
  while (currentdt >= dt)
    {
      currentdt -= dt;
      newG(PDat);
    }
}

void 
OPVACF::eventUpdate(const CIntEvent& iEvent, const C2ParticleData& PDat)
{
  //Move the time forward
  currentdt += iEvent.getdt();
  
  //Now test if we've gone over the step time
  while (currentdt >= dt)
    {
      currentdt -= dt;
      newG(PDat);
    }
}

void 
OPVACF::newG(const C1ParticleData& PDat)
{
  if (Sim->dynamics.liouvilleanTypeTest<LSLLOD>())
    Sim->dynamics.getLiouvillean().updateAllParticles();

  for (size_t i = 0; i < Sim->lN; ++i)
    G[i].push_front(Sim->vParticleList[i].getVelocity());	      
  
  //Now correct the fact that the wrong velocity has been pushed
  G[PDat.getParticle().getID()].front() = PDat.getOldVel();

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
OPVACF::newG(const C2ParticleData& PDat)
{
  for (size_t i = 0; i < Sim->lN; ++i)
    G[i].push_front(Sim->vParticleList[i].getVelocity());	      
  
  //Now correct the fact that the wrong velocity has been pushed
  G[PDat.particle1_.getParticle().getID()].front() 
    = PDat.particle1_.getOldVel();

  G[PDat.particle2_.getParticle().getID()].front() 
    = PDat.particle2_.getOldVel();

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
OPVACF::newG(const CNParticleData& PDat)
{  
  //This ensures the list stays at accumilator size
  for (size_t i = 0; i < Sim->lN; ++i)
    G[i].push_front(Sim->vParticleList[i].getVelocity());
  
  //Go back and fix the pushes
  BOOST_FOREACH(const C1ParticleData&PDat2, PDat.L1partChanges)
    G[PDat2.getParticle().getID()].front() = PDat2.getOldVel();		
  
  BOOST_FOREACH(const C2ParticleData& PDat2, PDat.L2partChanges)
    {
      G[PDat2.particle1_.getParticle().getID()].front() 
	= PDat2.particle1_.getOldVel();

      G[PDat2.particle2_.getParticle().getID()].front() 
	= PDat2.particle2_.getOldVel();
    }

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
OPVACF::output(xmlw::XmlStream& XML)
{
  Iflt factor = Sim->dynamics.units().unitTime() 
    / (Sim->dynamics.units().unitDiffusion() * count);

  for (size_t i = 0; i < accG2.size(); ++i)
    {
      Iflt specCount = Sim->dynamics.getSpecies()[i]->getCount();

      Vector  acc = 0.5*(accG2[i].front() + accG2[i].back());
      
      for (size_t j = 1; j < accG2[i].size() - 1; ++j)
	acc += accG2[i][j];
      
      acc *= factor * dt / (Sim->dynamics.units().unitTime() * specCount);

      XML << xmlw::tag("Correlator")
	  << xmlw::attr("name") << "VACF"
	  << xmlw::attr("species") << Sim->dynamics.getSpecies()[i]->getName()
	  << xmlw::attr("size") << accG2.size()
	  << xmlw::attr("dt") << dt / Sim->dynamics.units().unitTime()
	  << xmlw::attr("LengthInMFT") << dt * accG2[i].size() 
	/ Sim->getOutputPlugin<OPMisc>()->getMFT()
	  << xmlw::attr("simFactor") << factor / specCount
	  << xmlw::attr("SampleCount") << count
	  << xmlw::tag("Integral") << acc
	  << xmlw::endtag("Integral")
	  << xmlw::chardata();
            
      for (size_t j = 0; j < accG2[i].size(); ++j)
	{
	  XML << j * dt / Sim->dynamics.units().unitTime();
	  for (size_t iDim = 0; iDim < NDIM; iDim++)
	    XML << "\t" << accG2[i][j][iDim] * factor / specCount;
	  XML << "\n";
	}
      
      XML << xmlw::endtag("Correlator");
    }
}

Iflt 
OPVACF::getdt()
{
  //Get the simulation temperature
  if (dt == 0.0)
    {
      if (Sim->lastRunMFT != 0.0)
	return Sim->lastRunMFT * 50.0 / CorrelatorLength;
      else
	return 10.0 / (((Iflt) CorrelatorLength)*sqrt(Sim->dynamics.getLiouvillean().getkT()) * CorrelatorLength);
    }
  else 
    return dt;
}

void 
OPVACF::accPass()
{
  ++count;
  
  BOOST_FOREACH(const smrtPlugPtr<CSpecies>& spec, Sim->dynamics.getSpecies())
    BOOST_FOREACH(const size_t& ID, *spec->getRange())
    for (size_t j = 0; j < CorrelatorLength; ++j)
      for (size_t iDim(0); iDim < NDIM; ++iDim)
	accG2[spec->getID()][j][iDim] +=  G[ID].front()[iDim] * G[ID][j][iDim];
}
