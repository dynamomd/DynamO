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

#include "viscosityE.hpp"
#include <boost/foreach.hpp>
#include "../../dynamics/include.hpp"
#include "../0partproperty/misc.hpp"
#include "../1partproperty/kenergy.hpp"
#include "../../datatypes/vector.xml.hpp"

OPViscosityE::OPViscosityE(const DYNAMO::SimData* tmp, const XMLNode& XML):
  OutputPlugin(tmp,"ViscosityE", 60),
  count(0),
  dt(0),
  currentdt(0.0),
  currlen(0),
  notReady(true),
  CorrelatorLength(100),
  dtfactor(1.0)
{
  operator<<(XML);

  matrix zero;
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    for (size_t jDim = 0; jDim < NDIM; ++jDim)
      zero[iDim][jDim] = 0.0;

  G.resize(CorrelatorLength, zero);
    
  accG2.resize(CorrelatorLength, zero);

  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    for (size_t jDim = 0; jDim < NDIM; ++jDim)
      {
	avgTrace[iDim][jDim] = 0.0;
	delG[iDim][jDim] = 0.0;
	constDelG[iDim][jDim] = 0.0;
      }
}

void 
OPViscosityE::operator<<(const XMLNode& XML)
{
  try 
    {
      if (XML.isAttributeSet("Length"))
	CorrelatorLength = boost::lexical_cast<unsigned int>
	  (XML.getAttribute("Length"));

      if (XML.isAttributeSet("dt"))
	dt = Sim->dynamics.units().unitTime() * 
	  boost::lexical_cast<Iflt>(XML.getAttribute("dt"));

      if (XML.isAttributeSet("dtfactor"))
	dtfactor = boost::lexical_cast<Iflt>(XML.getAttribute("dtfactor"));
      
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
OPViscosityE::initialise()
{
  Sim->getOutputPlugin<OPKEnergy>();
  Sim->getOutputPlugin<OPMisc>();
  
  if (dt == 0.0)
    {
      if (Sim->lastRunMFT != 0.0)
	dt = Sim->lastRunMFT * 0.5 * dtfactor;
      else
	dt = 10.0 / (((Iflt) CorrelatorLength) * sqrt(Sim->dynamics.getLiouvillean().getkT()) * CorrelatorLength);
    }

  BOOST_FOREACH(const Particle& part, Sim->particleList)
    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      for (size_t jDim = 0; jDim < NDIM; ++jDim)
	constDelG[iDim][jDim] 
	  += part.getVelocity()[iDim] * part.getVelocity()[jDim]
	  * Sim->dynamics.getSpecies(part).getMass();

  I_cout() << "dt set to " << dt / Sim->dynamics.units().unitTime();
}

void 
OPViscosityE::eventUpdate(const GlobalEvent& iEvent, const NEventData& PDat) 
{
  stream(iEvent.getdt());
  impulseDelG(PDat);
  updateConstDelG(PDat);
}

void 
OPViscosityE::eventUpdate(const LocalEvent& iEvent, const NEventData& PDat) 
{
  stream(iEvent.getdt());
  impulseDelG(PDat);
  updateConstDelG(PDat);
}
  
void 
OPViscosityE::eventUpdate(const System&, const NEventData& PDat, const Iflt& edt) 
{ 
  stream(edt);
  impulseDelG(PDat);
  updateConstDelG(PDat);
}
  
void 
OPViscosityE::eventUpdate(const IntEvent& iEvent, const PairEventData& PDat)
{
  stream(iEvent.getdt());
  impulseDelG(PDat);
  updateConstDelG(PDat);
}

void 
OPViscosityE::stream(const Iflt& edt)
{
  //Move the time forward
  //currentdt += edt;
  
  //Now test if we've gone over the step time
  if (currentdt + edt >= dt)
    {
      for (size_t iDim = 0; iDim < NDIM; ++iDim)
	for (size_t jDim = 0; jDim < NDIM; ++jDim)
	  delG[iDim][jDim] += constDelG[iDim][jDim] * (dt - currentdt);

      newG (delG);

      currentdt += edt - dt;
      
      while (currentdt >= dt)
	{
	  for (size_t iDim = 0; iDim < NDIM; ++iDim)
	    for (size_t jDim = 0; jDim < NDIM; ++jDim)
	      delG[iDim][jDim] = constDelG[iDim][jDim] * dt;

	  currentdt -= dt;

	  newG(delG);
	}

      //Now calculate the start of the new delG
      for (size_t iDim = 0; iDim < NDIM; ++iDim)
	for (size_t jDim = 0; jDim < NDIM; ++jDim)
	  delG[iDim][jDim] = constDelG[iDim][jDim] * currentdt;

    }
  else
    {
      currentdt += edt;

      for (size_t iDim = 0; iDim < NDIM; ++iDim)
	for (size_t jDim = 0; jDim < NDIM; ++jDim)
	  delG[iDim][jDim] += constDelG[iDim][jDim] * edt;
    }
}

void 
OPViscosityE::newG(const matrix& Gval)
{
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    for (size_t jDim = 0; jDim < NDIM; ++jDim)
      avgTrace[iDim][jDim] += Gval[iDim][jDim];

  G.push_front(Gval);

  if (notReady)
    {
      if (++currlen != CorrelatorLength)
	return;

      notReady = false;
    }

  accPass();
}


void
OPViscosityE::impulseDelG(const PairEventData& colldat)
{
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    for (size_t jDim = 0; jDim < NDIM; ++jDim)
      delG[iDim][jDim] += colldat.particle1_.getDeltaP()[iDim] * colldat.rij[jDim];
}

void
OPViscosityE::impulseDelG(const NEventData& ndat) 
{ 
  BOOST_FOREACH(const PairEventData& dat, ndat.L2partChanges)
    for (size_t iDim(0); iDim < NDIM; ++iDim)
      for (size_t jDim(0); jDim < NDIM; ++jDim)
	delG[iDim][jDim] += dat.particle1_.getDeltaP()[iDim] 
	  * dat.rij[jDim];
}

inline void 
OPViscosityE::output(xml::XmlStream &XML)
{
  Iflt rescaleFactor = 1.0
    / (Sim->dynamics.units().unitTime() 
       //This line should be 1 however we have scaled the correlator time as well
       * Sim->dynamics.units().unitViscosity() * 2.0 
       * Sim->getOutputPlugin<OPKEnergy>()->getAvgkT() 
       //Count has been taken out due to the extra averaging of the constant piece 
       * Sim->dynamics.units().simVolume());
  
  XML << xml::tag("EinsteinCorrelator")
      << xml::attr("name") << name
      << xml::attr("size") << accG2.size()
      << xml::attr("dt") << dt / Sim->dynamics.units().unitTime()
      << xml::attr("LengthInMFT") << dt * accG2.size() / (Sim->getOutputPlugin<OPMisc>())->getMFT()
      << xml::attr("simFactor") << rescaleFactor
      << xml::attr("SampleCount") << count
      << xml::attr("columns")
      << "t ";
  
  char name[3] = "xx";
  
  for (size_t i = 0; i < NDIM; i++)
    for (size_t j = 0; j < NDIM; j++)
      {
	name[0] = 'x' + i;
	name[1] = 'x' + j;
	XML << name << " ";       
      }
  
  matrix traceAverage, P;

  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    for (size_t jDim = 0; jDim < NDIM; ++jDim)
      {
	traceAverage[iDim][jDim] = avgTrace[iDim][jDim] / (((Iflt) G.size()) + ((Iflt) count));
	
	P[iDim][jDim] = traceAverage[iDim][jDim] / (dt * Sim->dynamics.units().simVolume());
      }
  
  XML << xml::tag("Pressure");
 
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    {
      std::string name = std::string("d") + boost::lexical_cast<std::string>(iDim);
      
      XML << xml::tag(name.c_str());

      for (size_t jDim = 0; jDim < NDIM; ++jDim)
	{
	  std::string name = std::string("d") + boost::lexical_cast<std::string>(jDim);	  
	  XML << xml::attr(name.c_str())
	      << P[iDim][jDim] / Sim->dynamics.units().unitPressure();
	}
      
      XML << xml::endtag(name.c_str());
    }
  
  XML << xml::endtag("Pressure");
  
  Iflt AvgPressure = 0.0;
  for (size_t iDim = 0; iDim < NDIM; iDim++)
    AvgPressure += P[iDim][iDim];
  
  XML << xml::tag("PressureVals")
      << xml::attr("AvgPressure")
      << AvgPressure / (NDIM * Sim->dynamics.units().unitPressure())
      << xml::endtag("PressureVals");
  
  XML << xml::chardata();
  
  for (unsigned int i = 0; i < accG2.size(); i++)
    {
      XML << (i+1) * dt / Sim->dynamics.units().unitTime();
      for (size_t j = 0; j < NDIM; j++)
	for (size_t k = 0; k < NDIM; k++)
	  XML << "\t" << ((accG2[i][j][k] / count) - traceAverage[j][k]*(i+1)*traceAverage[j][k]*(i+1)) * rescaleFactor ;
      
      XML << "\n";
    }
  
  XML << xml::endtag("EinsteinCorrelator");
}

void 
OPViscosityE::updateConstDelG(const PairEventData& PDat)
{
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    for (size_t jDim = 0; jDim < NDIM; ++jDim)
      constDelG[iDim][jDim] += 
	(PDat.particle1_.getParticle().getVelocity()[iDim]
	 * PDat.particle1_.getParticle().getVelocity()[jDim]
	 - PDat.particle1_.getOldVel()[iDim]
	 * PDat.particle1_.getOldVel()[jDim])
	* PDat.particle1_.getSpecies().getMass()
	+ (PDat.particle2_.getParticle().getVelocity()[iDim]
	   * PDat.particle2_.getParticle().getVelocity()[jDim]
	   - PDat.particle2_.getOldVel()[iDim]
	   * PDat.particle2_.getOldVel()[jDim])
	* PDat.particle2_.getSpecies().getMass();
}

void 
OPViscosityE::updateConstDelG(const ParticleEventData& PDat)
{
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    for (size_t jDim = 0; jDim < NDIM; ++jDim)
      constDelG[iDim][jDim] += 
	(PDat.getParticle().getVelocity()[iDim]
	 * PDat.getParticle().getVelocity()[jDim]
	 - PDat.getOldVel()[iDim] * PDat.getOldVel()[jDim]
	 ) * PDat.getSpecies().getMass();
}

void 
OPViscosityE::updateConstDelG(const NEventData& ndat)
{
  BOOST_FOREACH(const ParticleEventData& dat, ndat.L1partChanges)
    updateConstDelG(dat);
  
  BOOST_FOREACH(const PairEventData& dat, ndat.L2partChanges)
    updateConstDelG(dat);
}

void 
OPViscosityE::accPass()
{
  ++count;
  matrix sum;

  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    for (size_t jDim = 0; jDim < NDIM; ++jDim)
      sum[iDim][jDim] = 0.0;

  for (size_t i = 0; i < CorrelatorLength; ++i)
    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      for (size_t jDim = 0; jDim < NDIM; ++jDim)
	{
	  sum[iDim][jDim] += G[i][iDim][jDim];
	  accG2[i][iDim][jDim] += sum[iDim][jDim] * sum[iDim][jDim];
	}
}
