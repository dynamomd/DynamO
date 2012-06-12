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

#include <dynamo/outputplugins/correlations/viscosityE.hpp>
#include <dynamo/dynamics/include.hpp>
#include <dynamo/outputplugins/0partproperty/misc.hpp>
#include <dynamo/outputplugins/1partproperty/kenergy.hpp>
#include <boost/foreach.hpp>

namespace dynamo {
  OPViscosityE::OPViscosityE(const dynamo::SimData* tmp, 
			     const magnet::xml::Node& XML):
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
  OPViscosityE::operator<<(const magnet::xml::Node& XML)
  {
    try 
      {
	//The length of the correlation function in timesteps
	if (XML.hasAttribute("Length"))
	  CorrelatorLength = XML.getAttribute("Length").as<size_t>();

	//The time step between sampling the correlation function
	if (XML.hasAttribute("dt"))
	  dt = Sim->dynamics.units().unitTime() * 
	    XML.getAttribute("dt").as<double>();

	//This is used to set the timestep as a function of the previous
	//simulation runs mean free time.
	if (XML.hasAttribute("dtfactor"))
	  dtfactor = XML.getAttribute("dtfactor").as<double>();
      
	//This sets the total correlation time, and the time step is
	//worked out as \f$ dt = t / CorrelatorLength \f$.
	if (XML.hasAttribute("t"))
	  dt = Sim->dynamics.units().unitTime() * 
	    XML.getAttribute("t").as<double>() / CorrelatorLength;
      }
    catch (boost::bad_lexical_cast &)
      {
	M_throw() << "Failed a lexical cast in OPVACF";
      }  
  }

  void 
  OPViscosityE::initialise()
  {
    if (!(Sim->getOutputPlugin<OPMisc>()))
      M_throw() << "ViscosityE requires Misc output plugin!";
    if (!(Sim->getOutputPlugin<OPKEnergy>()))
      M_throw() << "ViscosityE requires KEnergy output plugin!";
  
    if (dt == 0.0)
      {
	if (Sim->lastRunMFT != 0.0)
	  dt = Sim->lastRunMFT * 0.5 * dtfactor;
	else
	  dt = 10.0 / (((double) CorrelatorLength) * sqrt(Sim->dynamics.getLiouvillean().getkT()) * CorrelatorLength);
      }

    BOOST_FOREACH(const Particle& part, Sim->particleList)
      for (size_t iDim = 0; iDim < NDIM; ++iDim)
	for (size_t jDim = 0; jDim < NDIM; ++jDim)
	  constDelG[iDim][jDim] 
	    += part.getVelocity()[iDim] * part.getVelocity()[jDim]
	    * Sim->species[part].getMass(part.getID());

    dout << "dt set to " << dt / Sim->dynamics.units().unitTime() << std::endl;
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
  OPViscosityE::eventUpdate(const System&, const NEventData& PDat, const double& edt) 
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
  OPViscosityE::stream(const double& edt)
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
  OPViscosityE::output(magnet::xml::XmlStream &XML)
  {
    double rescaleFactor = 1.0
      / (Sim->dynamics.units().unitTime() 
	 //This line should be 1 however we have scaled the correlator time as well
	 * Sim->dynamics.units().unitViscosity() * 2.0 
	 * Sim->getOutputPlugin<OPKEnergy>()->getAvgkT() 
	 //Count has been taken out due to the extra averaging of the constant piece 
	 * Sim->dynamics.getSimVolume());
  
    XML << magnet::xml::tag("EinsteinCorrelator")
	<< magnet::xml::attr("name") << name
	<< magnet::xml::attr("size") << accG2.size()
	<< magnet::xml::attr("dt") << dt / Sim->dynamics.units().unitTime()
	<< magnet::xml::attr("LengthInMFT") << dt * accG2.size() / (Sim->getOutputPlugin<OPMisc>())->getMFT()
	<< magnet::xml::attr("simFactor") << rescaleFactor
	<< magnet::xml::attr("SampleCount") << count
	<< magnet::xml::attr("columns")
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
	
	  P[iDim][jDim] = traceAverage[iDim][jDim] / (dt * Sim->dynamics.getSimVolume());
	}
  
    XML << magnet::xml::tag("Pressure");
 
    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      {
	std::string name = std::string("d") + boost::lexical_cast<std::string>(iDim);
      
	XML << magnet::xml::tag(name.c_str());

	for (size_t jDim = 0; jDim < NDIM; ++jDim)
	  {
	    std::string name = std::string("d") + boost::lexical_cast<std::string>(jDim);	  
	    XML << magnet::xml::attr(name.c_str())
		<< P[iDim][jDim] / Sim->dynamics.units().unitPressure();
	  }
      
	XML << magnet::xml::endtag(name.c_str());
      }
  
    XML << magnet::xml::endtag("Pressure");
  
    double AvgPressure = 0.0;
    for (size_t iDim = 0; iDim < NDIM; iDim++)
      AvgPressure += P[iDim][iDim];
  
    XML << magnet::xml::tag("PressureVals")
	<< magnet::xml::attr("AvgPressure")
	<< AvgPressure / (NDIM * Sim->dynamics.units().unitPressure())
	<< magnet::xml::endtag("PressureVals");
  
    XML << magnet::xml::chardata();
  
    for (unsigned int i = 0; i < accG2.size(); i++)
      {
	XML << (i+1) * dt / Sim->dynamics.units().unitTime();
	for (size_t j = 0; j < NDIM; j++)
	  for (size_t k = 0; k < NDIM; k++)
	    XML << "\t" << ((accG2[i][j][k] / count) - traceAverage[j][k]*(i+1)*traceAverage[j][k]*(i+1)) * rescaleFactor ;
      
	XML << "\n";
      }
  
    XML << magnet::xml::endtag("EinsteinCorrelator");
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
	  * PDat.particle1_.getSpecies().getMass(PDat.particle1_.getParticle().getID())
	  + (PDat.particle2_.getParticle().getVelocity()[iDim]
	     * PDat.particle2_.getParticle().getVelocity()[jDim]
	     - PDat.particle2_.getOldVel()[iDim]
	     * PDat.particle2_.getOldVel()[jDim])
	  * PDat.particle2_.getSpecies().getMass(PDat.particle2_.getParticle().getID());
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
	   ) * PDat.getSpecies().getMass(PDat.getParticle().getID());
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
}
