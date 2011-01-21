/*  DYNAMO:- Event driven molecular dynamics simulator 
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

#include "viscosityCollisionalE.hpp"
#include <boost/foreach.hpp>
#include "../../dynamics/include.hpp"
#include "../0partproperty/misc.hpp"
#include "../1partproperty/kenergy.hpp"
#include "../../datatypes/vector.xml.hpp"

OPViscosityCollisionalE::OPViscosityCollisionalE(const DYNAMO::SimData* tmp, const XMLNode& XML):
  OutputPlugin(tmp,"ViscosityCollisionalE", 60),
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
OPViscosityCollisionalE::operator<<(const XMLNode& XML)
{
  try 
    {
      if (XML.isAttributeSet("Length"))
	CorrelatorLength = boost::lexical_cast<unsigned int>
	  (XML.getAttribute("Length"));

      if (XML.isAttributeSet("dt"))
	dt = Sim->dynamics.units().unitTime() * 
	  boost::lexical_cast<double>(XML.getAttribute("dt"));

      if (XML.isAttributeSet("dtfactor"))
	dtfactor = boost::lexical_cast<double>(XML.getAttribute("dtfactor"));
      
      if (XML.isAttributeSet("t"))
	dt = Sim->dynamics.units().unitTime() * 
	  boost::lexical_cast<double>(XML.getAttribute("t"))/CorrelatorLength;
    }
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in OPVACF";
    }  
}

void 
OPViscosityCollisionalE::initialise()
{
  Sim->getOutputPlugin<OPMisc>();
  
  if (dt == 0.0)
    {
      if (Sim->lastRunMFT != 0.0)
	dt = Sim->lastRunMFT * 0.5 * dtfactor;
      else
	dt = 10.0 / (((double) CorrelatorLength) * sqrt(Sim->dynamics.getLiouvillean().getkT()) * CorrelatorLength);
    }

  I_cout() << "dt set to " << dt / Sim->dynamics.units().unitTime();
}

void 
OPViscosityCollisionalE::eventUpdate(const GlobalEvent& iEvent, 
				      const NEventData& PDat) 
{
  stream(iEvent.getdt());
  impulseDelG(PDat);
}

void 
OPViscosityCollisionalE::eventUpdate(const LocalEvent& iEvent, 
				      const NEventData& PDat) 
{
  stream(iEvent.getdt());
  impulseDelG(PDat);
}
  
void 
OPViscosityCollisionalE::eventUpdate(const System&, 
				      const NEventData& PDat, 
				      const double& edt) 
{ 
  stream(edt);
  impulseDelG(PDat);
}
  
void 
OPViscosityCollisionalE::eventUpdate(const IntEvent& iEvent, 
				      const PairEventData& PDat)
{
  stream(iEvent.getdt());
  impulseDelG(PDat);
}

void 
OPViscosityCollisionalE::stream(const double& edt)
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
OPViscosityCollisionalE::newG(const matrix& Gval)
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
OPViscosityCollisionalE::impulseDelG(const PairEventData& colldat)
{
  for (size_t iDim = 0; iDim < NDIM; ++iDim)
    for (size_t jDim = 0; jDim < NDIM; ++jDim)
      delG[iDim][jDim] += colldat.particle1_.getDeltaP()[iDim] * colldat.rij[jDim];
}

void
OPViscosityCollisionalE::impulseDelG(const NEventData& ndat) 
{ 
  BOOST_FOREACH(const PairEventData& dat, ndat.L2partChanges)
    for (size_t iDim(0); iDim < NDIM; ++iDim)
      for (size_t jDim(0); jDim < NDIM; ++jDim)
	delG[iDim][jDim] += dat.particle1_.getDeltaP()[iDim] 
	  * dat.rij[jDim];
}

inline void 
OPViscosityCollisionalE::output(xml::XmlStream &XML)
{
  double rescaleFactor = 1.0
    / (Sim->dynamics.units().unitTime() 
       //This line should be 1 however we have scaled the correlator time as well
       * Sim->dynamics.units().unitViscosity() * 2.0 
       //Count has been taken out due to the extra averaging of the constant piece 
       * Sim->dynamics.units().simVolume());
  
  XML << xml::tag("EinsteinCorrelator")
      << xml::attr("name") << "ViscosityTimesT"
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
	traceAverage[iDim][jDim] = avgTrace[iDim][jDim] / (((double) G.size()) + ((double) count));
	
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
  
  double AvgPressure = 0.0;
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
	  if (k==j)
	    XML << "\t" << ((accG2[i][j][k] / count) - pow(traceAverage[j][k]*(i+1),2)) * rescaleFactor ;
	  else
	    XML << "\t" << ((accG2[i][j][k] / count)) * rescaleFactor; 
      
      XML << "\n";
    }
  
  XML << xml::endtag("EinsteinCorrelator");
}

void 
OPViscosityCollisionalE::accPass()
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
