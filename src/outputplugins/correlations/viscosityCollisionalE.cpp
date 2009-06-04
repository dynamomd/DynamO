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

#include "viscosityCollisionalE.hpp"
#include <boost/foreach.hpp>
#include "../../dynamics/include.hpp"
#include "../0partproperty/misc.hpp"
#include "../1partproperty/kenergy.hpp"
#include "../../datatypes/vector.xml.hpp"

COPViscosityCollisionalE::COPViscosityCollisionalE(const DYNAMO::SimData* tmp, const XMLNode& XML):
  COutputPlugin(tmp,"ViscosityCollisionalE", 60),
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
      }
}

void 
COPViscosityCollisionalE::operator<<(const XMLNode& XML)
{
  try 
    {
      if (XML.isAttributeSet("Length"))
	CorrelatorLength = boost::lexical_cast<unsigned int>
	  (XML.getAttribute("Length"));

      if (XML.isAttributeSet("dt"))
	dt = Sim->Dynamics.units().unitTime() * 
	  boost::lexical_cast<Iflt>(XML.getAttribute("dt"));

      if (XML.isAttributeSet("dtfactor"))
	dtfactor = boost::lexical_cast<Iflt>(XML.getAttribute("dtfactor"));
      
      if (XML.isAttributeSet("t"))
	dt = Sim->Dynamics.units().unitTime() * 
	  boost::lexical_cast<Iflt>(XML.getAttribute("t"))/CorrelatorLength;
    }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in COPVACF";
    }  
}

void 
COPViscosityCollisionalE::initialise()
{
  Sim->getOutputPlugin<COPMisc>();
  
  if (dt == 0.0)
    {
      if (Sim->lastRunMFT != 0.0)
	dt = Sim->lastRunMFT * 0.5 * dtfactor;
      else
	dt = 10.0 / (((Iflt) CorrelatorLength) * sqrt(Sim->Dynamics.Liouvillean().getkT()) * CorrelatorLength);
    }

  I_cout() << "dt set to " << dt / Sim->Dynamics.units().unitTime();
}

void 
COPViscosityCollisionalE::eventUpdate(const CGlobEvent& iEvent, 
				      const CNParticleData& PDat) 
{
  stream(iEvent.getdt());
  impulseDelG(PDat);
}

void 
COPViscosityCollisionalE::eventUpdate(const CLocalEvent& iEvent, 
				      const CNParticleData& PDat) 
{
  stream(iEvent.getdt());
  impulseDelG(PDat);
}
  
void 
COPViscosityCollisionalE::eventUpdate(const CSystem&, 
				      const CNParticleData& PDat, 
				      const Iflt& edt) 
{ 
  stream(edt);
  impulseDelG(PDat);
}
  
void 
COPViscosityCollisionalE::eventUpdate(const CIntEvent& iEvent, 
				      const C2ParticleData& PDat)
{
  stream(iEvent.getdt());
  impulseDelG(PDat);
}

void 
COPViscosityCollisionalE::stream(const Iflt& edt)
{
  //Move the time forward
  //currentdt += edt;
  
  //Now test if we've gone over the step time
  if (currentdt + edt >= dt)
    {
      newG (delG);

      currentdt += edt - dt;

      //Now calculate the start of the new delG
      for (size_t iDim = 0; iDim < NDIM; ++iDim)
	for (size_t jDim = 0; jDim < NDIM; ++jDim)
	  delG[iDim][jDim] = 0.0;
      
      while (currentdt >= dt)
	{
	  currentdt -= dt;

	  newG(delG);
	}
    }
  else
      currentdt += edt;
}

void 
COPViscosityCollisionalE::newG(const matrix& Gval)
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
COPViscosityCollisionalE::impulseDelG(const C2ParticleData& colldat)
{
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    for (size_t jDim = 0; jDim < NDIM; ++jDim)
      delG[iDim][jDim] += colldat.particle1_.getDeltaP()[iDim] * colldat.rij[jDim];
}

void
COPViscosityCollisionalE::impulseDelG(const CNParticleData& ndat) 
{ 
  BOOST_FOREACH(const C2ParticleData& dat, ndat.L2partChanges)
    for (size_t iDim(0); iDim < NDIM; ++iDim)
      for (size_t jDim(0); jDim < NDIM; ++jDim)
	delG[iDim][jDim] += dat.particle1_.getDeltaP()[iDim] 
	  * dat.rij[jDim];
}

inline void 
COPViscosityCollisionalE::output(xmlw::XmlStream &XML)
{
  Iflt rescaleFactor = 1.0
    / (Sim->Dynamics.units().unitTime() 
       //This line should be 1 however we have scaled the correlator time as well
       * Sim->Dynamics.units().unitViscosity() * 2.0 
       //Count has been taken out due to the extra averaging of the constant piece 
       * Sim->Dynamics.units().simVolume());
  
  XML << xmlw::tag("EinsteinCorrelator")
      << xmlw::attr("name") << "ViscosityTimesT"
      << xmlw::attr("size") << accG2.size()
      << xmlw::attr("dt") << dt / Sim->Dynamics.units().unitTime()
      << xmlw::attr("LengthInMFT") << dt * accG2.size() / (Sim->getOutputPlugin<COPMisc>())->getMFT()
      << xmlw::attr("simFactor") << rescaleFactor
      << xmlw::attr("SampleCount") << count
      << xmlw::attr("columns")
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
	
	P[iDim][jDim] = traceAverage[iDim][jDim] / (dt * Sim->Dynamics.units().simVolume());
      }
  
  XML << xmlw::tag("Pressure");
 
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    {
      std::string name = std::string("d") + boost::lexical_cast<std::string>(iDim);
      
      XML << xmlw::tag(name.c_str());

      for (size_t jDim = 0; jDim < NDIM; ++jDim)
	{
	  std::string name = std::string("d") + boost::lexical_cast<std::string>(jDim);	  
	  XML << xmlw::attr(name.c_str())
	      << P[iDim][jDim] / Sim->Dynamics.units().unitPressure();
	}
      
      XML << xmlw::endtag(name.c_str());
    }
  
  XML << xmlw::endtag("Pressure");
  
  Iflt AvgPressure = 0.0;
  for (size_t iDim = 0; iDim < NDIM; iDim++)
    AvgPressure += P[iDim][iDim];
  
  XML << xmlw::tag("PressureVals")
      << xmlw::attr("AvgPressure")
      << AvgPressure / (NDIM * Sim->Dynamics.units().unitPressure())
      << xmlw::endtag("PressureVals");
  
  XML << xmlw::chardata();
  
  for (unsigned int i = 0; i < accG2.size(); i++)
    {
      XML << (i+1) * dt / Sim->Dynamics.units().unitTime();
      for (size_t j = 0; j < NDIM; j++)
	for (size_t k = 0; k < NDIM; k++)
	  if (k==j)
	    XML << "\t" << ((accG2[i][j][k] / count) - pow(traceAverage[j][k]*(i+1),2)) * rescaleFactor ;
	  else
	    XML << "\t" << ((accG2[i][j][k] / count)) * rescaleFactor; 
      
      XML << "\n";
    }
  
  XML << xmlw::endtag("EinsteinCorrelator");
}

void 
COPViscosityCollisionalE::accPass()
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
