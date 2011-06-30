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

#include "RingDSMC.hpp"
#include "../dynamics.hpp"
#include "../units/units.hpp"
#include "../BC/BC.hpp"
#include "../../simulation/particle.hpp"
#include "../species/species.hpp"
#include "../NparticleEventData.hpp"
#include "../ranges/include.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../../schedulers/scheduler.hpp"
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>
#include <boost/random/uniform_int.hpp>

#ifdef DYNAMO_DEBUG 
#include <boost/math/special_functions/fpclassify.hpp>
#endif

CSRingDSMC::CSRingDSMC(const magnet::xml::Node& XML, dynamo::SimData* tmp): 
  System(tmp),
  uniformRand(Sim->ranGenerator, boost::uniform_real<>(0,1)),
  maxprob12(0.0),
  maxprob13(0.0),
  range1(NULL)
{
  dt = HUGE_VAL;
  operator<<(XML);
  type = DSMC;
}

CSRingDSMC::CSRingDSMC(dynamo::SimData* nSim, double nd, double ntstp, double nChi1, double nChi2,
			     double ne, std::string nName, CRange* r1):
  System(nSim),
  uniformRand(Sim->ranGenerator,boost::uniform_real<>(0,1)),
  tstep(ntstp),
  chi12(nChi1),
  chi13(nChi2),
  d2(nd * nd),
  diameter(nd),
  maxprob12(0.0),
  maxprob13(0.0),
  e(ne),
  range1(r1)
{
  sysName = nName;
  type = DSMC;
  if (r1->size() % 2)
    M_throw() << "Need an even number of particles in"
	      << " the range to make a whole number of velocity pairs";
}

void 
CSRingDSMC::runEvent() const
{
  double locdt = dt;
  
#ifdef DYNAMO_DEBUG 
  if (boost::math::isnan(locdt))
    M_throw() << "A NAN system event time has been found";
#endif
    
  Sim->dSysTime += locdt;
    
  Sim->ptrScheduler->stream(locdt);
  
  //dynamics must be updated first
  Sim->dynamics.stream(locdt);

  dt = tstep;

  locdt += Sim->freestreamAcc;
  Sim->freestreamAcc = 0;

  BOOST_FOREACH(magnet::ClonePtr<OutputPlugin>& Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(*this, NEventData(), locdt);

  //////////////////// T(1,2) operator
  double intPart;
  double fracpart = std::modf(maxprob12 * range1->size(), &intPart);
 
  size_t nmax = static_cast<size_t>(intPart) + (Sim->uniform_sampler() < fracpart);
  
  {
    boost::variate_generator
      <dynamo::baseRNG&, boost::uniform_int<size_t> >
      id1sampler(Sim->ranGenerator, 
		 boost::uniform_int<size_t>(0, (range1->size()/2) - 1));
    
    for (size_t n = 0; n < nmax; ++n)
      {
	size_t pairID(id1sampler());
	const Particle& p1(Sim->particleList[*(range1->begin() + 2 * pairID)]);
	const Particle& p2(Sim->particleList[*(range1->begin() + 2 * pairID + 1)]);
	
	Sim->dynamics.getLiouvillean().updateParticlePair(p1, p2);
	
	CPDData PDat;
	
	for (size_t iDim(0); iDim < NDIM; ++iDim)
	  PDat.rij[iDim] = Sim->normal_sampler();
	
	PDat.rij *= diameter / PDat.rij.nrm();
	
	if (Sim->dynamics.getLiouvillean().DSMCSpheresTest
	    (p1, p2, maxprob12, factor12, PDat))
	  {
	    ++Sim->eventCount;
	    ++n12;

	    const PairEventData
	      SDat(Sim->dynamics.getLiouvillean().DSMCSpheresRun(p1, p2, e, PDat));
	    
	    Sim->signalParticleUpdate(SDat);
	    
	    Sim->ptrScheduler->fullUpdate(p1, p2);
	    
	    BOOST_FOREACH(magnet::ClonePtr<OutputPlugin>& Ptr, Sim->outputPlugins)
	      Ptr->eventUpdate(*this, SDat, 0.0);
	  }
      }
  }

  //////////////////// T(1,3) operator
  {
    fracpart = std::modf(maxprob13 * range1->size(), &intPart);
    
    nmax = static_cast<size_t>(intPart) + (Sim->uniform_sampler() < fracpart);
    
    boost::variate_generator
      <dynamo::baseRNG&, boost::uniform_int<size_t> >
      id1sampler(Sim->ranGenerator, 
		 boost::uniform_int<size_t>(0, range1->size() - 1));
    
    for (size_t n = 0; n < nmax; ++n)
      {
	const Particle& p1(Sim->particleList[*(range1->begin() + id1sampler())]);
	
	size_t secondID(id1sampler());
	
	while ((secondID == p1.getID())
	       || ((secondID % 2) 
		   ? ((secondID-1) == p1.getID())
		   : ((secondID+1) == p1.getID())))
	  secondID = id1sampler();
	
	const Particle& p2(Sim->particleList[*(range1->begin() + secondID)]);
	
	Sim->dynamics.getLiouvillean().updateParticlePair(p1, p2);
	
	CPDData PDat;
	
	for (size_t iDim(0); iDim < NDIM; ++iDim)
	  PDat.rij[iDim] = Sim->normal_sampler();
	
	PDat.rij *= diameter / PDat.rij.nrm();
	
	if (Sim->dynamics.getLiouvillean().DSMCSpheresTest
	    (p1, p2, maxprob13, factor13, PDat))
	  {
	    ++Sim->eventCount;
	    ++n13;

	    const PairEventData
	      SDat(Sim->dynamics.getLiouvillean().DSMCSpheresRun(p1, p2, e, PDat));

	    Sim->signalParticleUpdate(SDat);
	    
	    Sim->ptrScheduler->fullUpdate(p1, p2);
	    
	    BOOST_FOREACH(magnet::ClonePtr<OutputPlugin>& Ptr, Sim->outputPlugins)
	      Ptr->eventUpdate(*this, SDat, 0.0);
	  }
      }
  }
}

void
CSRingDSMC::initialise(size_t nID)
{
  ID = nID;
  dt = tstep;

  n12 = 0;
  n13 = 0;

  factor12 = range1->size()
    * diameter * M_PI * chi12 * tstep 
    / Sim->dynamics.getSimVolume();
  
  factor13 = range1->size()
    * diameter * M_PI * chi13 * tstep 
    / Sim->dynamics.getSimVolume();
  
  if (maxprob12 == 0.0)
    { 
      boost::variate_generator
	<dynamo::baseRNG&, boost::uniform_int<size_t> >
	id1sampler(Sim->ranGenerator, 
		   boost::uniform_int<size_t>(0, (range1->size()/2) - 1));
      
      //Just do some quick testing to get an estimate
      for (size_t n = 0; n < 1000; ++n)
	{
	  size_t pairID(id1sampler());
	  const Particle& p1(Sim->particleList[*(range1->begin() + 2 * pairID)]);
	  const Particle& p2(Sim->particleList[*(range1->begin() + 2 * pairID + 1)]);
	  
	  Sim->dynamics.getLiouvillean().updateParticlePair(p1, p2);
	  
	  CPDData PDat;
	  
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    PDat.rij[iDim] = Sim->normal_sampler();
	  
	  PDat.rij *= diameter / PDat.rij.nrm();
	  
	  Sim->dynamics.getLiouvillean().DSMCSpheresTest(p1, p2, maxprob12, 
						      factor12, PDat);
	}
    }

  if (maxprob13 == 0.0)
    { 
      boost::variate_generator
	<dynamo::baseRNG&, boost::uniform_int<size_t> >
	id1sampler(Sim->ranGenerator, 
		   boost::uniform_int<size_t>(0, range1->size() - 1));
      
      //Just do some quick testing to get an estimate
      for (size_t n = 0; n < 1000; ++n)
	{
	  const Particle& p1(Sim->particleList[*(range1->begin() + id1sampler())]);

	  size_t secondID(id1sampler());

	  while ((secondID == p1.getID())
		 || ((secondID % 2) 
		     ? ((secondID-1) == p1.getID())
		     : ((secondID+1) == p1.getID())))
	    secondID = id1sampler();

	  const Particle& p2(Sim->particleList[*(range1->begin() + secondID)]);
	  
	  Sim->dynamics.getLiouvillean().updateParticlePair(p1, p2);
	  
	  CPDData PDat;
	  
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    PDat.rij[iDim] = Sim->normal_sampler();
	  
	  PDat.rij *= diameter / PDat.rij.nrm();
	  
	  Sim->dynamics.getLiouvillean().DSMCSpheresTest(p1, p2, maxprob13, 
						      factor13, PDat);
	}
    }

  if (maxprob12 > 0.5)
    derr << "MaxProbability12 is " << maxprob12
	     << "\nNpairs12 per step is " << range1->size() * maxprob12 << std::endl;
  else
    dout << "MaxProbability12 is " << maxprob12
	     << "\nNpairs12 per step is " << range1->size() * maxprob12 << std::endl;

  if (maxprob13 > 0.5)
    derr << "MaxProbability13 is " << maxprob13
	     << "\nNpairs13 per step is " << range1->size() * maxprob13 << std::endl;
  else
    dout << "MaxProbability13 is " << maxprob13
	     << "\nNpairs13 per step is " << range1->size() * maxprob13 << std::endl;
  
  if (range1->size() * maxprob12 < 2.0)
    derr << "The 12 probability is low" << std::endl;

  if (range1->size() * maxprob13 < 2.0)
    derr << "The 13 probability is low" << std::endl;
}

void
CSRingDSMC::operator<<(const magnet::xml::Node& XML)
{
  if (strcmp(XML.getAttribute("Type"),"RingDSMC"))
    M_throw() << "Attempting to load RingDSMC from a " 
	      << XML.getAttribute("Type") <<  " entry"; 
  
  try {
    tstep = XML.getAttribute("tStep").as<double>() * Sim->dynamics.units().unitTime();    
    chi12 = XML.getAttribute("Chi12").as<double>();
    chi13 = XML.getAttribute("Chi13").as<double>();
    sysName = XML.getAttribute("Name");
    diameter = XML.getAttribute("Diameter").as<double>() * Sim->dynamics.units().unitLength();
    e = XML.getAttribute("Inelasticity").as<double>();
    d2 = diameter * diameter;
    range1.set_ptr(CRange::getClass(XML.getNode("Range1"), Sim));

    if (XML.hasAttribute("MaxProbability12"))
      maxprob12 = XML.getAttribute("MaxProbability12").as<double>();
	
    if (XML.hasAttribute("MaxProbability13"))
      maxprob13 = XML.getAttribute("MaxProbability13").as<double>();
  }
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in CGGlobal";
    }
}

void 
CSRingDSMC::outputXML(magnet::xml::XmlStream& XML) const
{
  if (n12 || n13)
    dout << "Number of T(1,2) events " << n12
	 << "\nNumber of T(1,3) events " << n13
	 << "\nRatio T(1,2)/total " << ((double) n12) / (((double) n13) + ((double) n12))
	 << std::endl;

  XML << magnet::xml::tag("System")
      << magnet::xml::attr("Type") << "RingDSMC"
      << magnet::xml::attr("tStep") << tstep / Sim->dynamics.units().unitTime()
      << magnet::xml::attr("Chi12") << chi12
      << magnet::xml::attr("Chi13") << chi13
      << magnet::xml::attr("Diameter") << diameter / Sim->dynamics.units().unitLength()
      << magnet::xml::attr("Inelasticity") << e
      << magnet::xml::attr("Name") << sysName
      << magnet::xml::attr("MaxProbability12") << maxprob12
      << magnet::xml::attr("MaxProbability13") << maxprob13
      << magnet::xml::tag("Range1")
      << range1
      << magnet::xml::endtag("Range1")
      << magnet::xml::endtag("System");
}
