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

#include <dynamo/systems/RingDSMC.hpp>

#include <dynamo/units/units.hpp>
#include <dynamo/BC/BC.hpp>
#include <dynamo/particle.hpp>
#include <dynamo/species/species.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/ranges/include.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>
#include <boost/random/uniform_int.hpp>

#ifdef DYNAMO_DEBUG 
#include <boost/math/special_functions/fpclassify.hpp>
#endif

namespace dynamo {
  SysRingDSMC::SysRingDSMC(const magnet::xml::Node& XML, dynamo::Simulation* tmp): 
    System(tmp),
    uniformRand(Sim->ranGenerator, boost::uniform_real<>(0,1)),
    maxprob12(0.0),
    maxprob13(0.0)
  {
    dt = HUGE_VAL;
    operator<<(XML);
    type = DSMC;
  }

  SysRingDSMC::SysRingDSMC(dynamo::Simulation* nSim, double nd, double ntstp, double nChi1, double nChi2,
			 double ne, std::string nName, Range* r1):
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
  SysRingDSMC::runEvent() const
  {
    double locdt = dt;
  
#ifdef DYNAMO_DEBUG 
    if (boost::math::isnan(locdt))
      M_throw() << "A NAN system event time has been found";
#endif
    
    Sim->systemTime += locdt;
    
    Sim->ptrScheduler->stream(locdt);
  
    //dynamics must be updated first
    Sim->stream(locdt);

    dt = tstep;

    BOOST_FOREACH(shared_ptr<OutputPlugin>& Ptr, Sim->outputPlugins)
      Ptr->eventUpdate(*this, NEventData(), locdt);

    //////////////////// T(1,2) operator
    double Event;
    double fracpart = std::modf(maxprob12 * range1->size(), &Event);
 
    size_t nmax = static_cast<size_t>(Event) + (Sim->uniform_sampler() < fracpart);
  
    {
      boost::variate_generator
	<dynamo::baseRNG&, boost::uniform_int<size_t> >
	id1sampler(Sim->ranGenerator, 
		   boost::uniform_int<size_t>(0, (range1->size()/2) - 1));
    
      for (size_t n = 0; n < nmax; ++n)
	{
	  size_t pairID(id1sampler());
	  Particle& p1(Sim->particles[*(range1->begin() + 2 * pairID)]);
	  Particle& p2(Sim->particles[*(range1->begin() + 2 * pairID + 1)]);
	
	  Sim->dynamics->updateParticlePair(p1, p2);
	
	  Vector rij;
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    rij[iDim] = Sim->normal_sampler();
	  rij *= diameter / rij.nrm();
	
	  if (Sim->dynamics->DSMCSpheresTest
	      (p1, p2, maxprob12, factor12, rij))
	    {
	      ++Sim->eventCount;
	      ++n12;

	      const PairEventData 
		SDat(Sim->dynamics->DSMCSpheresRun(p1, p2, e, rij));
	    
	      Sim->signalParticleUpdate(SDat);
	    
	      Sim->ptrScheduler->fullUpdate(p1, p2);
	    
	      BOOST_FOREACH(shared_ptr<OutputPlugin>& Ptr, Sim->outputPlugins)
		Ptr->eventUpdate(*this, SDat, 0.0);
	    }
	}
    }

    //////////////////// T(1,3) operator
    {
      fracpart = std::modf(maxprob13 * range1->size(), &Event);
    
      nmax = static_cast<size_t>(Event) + (Sim->uniform_sampler() < fracpart);
    
      boost::variate_generator
	<dynamo::baseRNG&, boost::uniform_int<size_t> >
	id1sampler(Sim->ranGenerator, 
		   boost::uniform_int<size_t>(0, range1->size() - 1));
    
      for (size_t n = 0; n < nmax; ++n)
	{
	  Particle& p1(Sim->particles[*(range1->begin() + id1sampler())]);
	
	  size_t secondID(id1sampler());
	
	  while ((secondID == p1.getID())
		 || ((secondID % 2) 
		     ? ((secondID-1) == p1.getID())
		     : ((secondID+1) == p1.getID())))
	    secondID = id1sampler();
	
	  Particle& p2(Sim->particles[*(range1->begin() + secondID)]);
	
	  Sim->dynamics->updateParticlePair(p1, p2);
	
	  Vector rij;
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    rij[iDim] = Sim->normal_sampler();
	
	  rij *= diameter / rij.nrm();
	
	  if (Sim->dynamics->DSMCSpheresTest
	      (p1, p2, maxprob13, factor13, rij))
	    {
	      ++Sim->eventCount;
	      ++n13;

	      const PairEventData
		SDat(Sim->dynamics->DSMCSpheresRun(p1, p2, e, rij));

	      Sim->signalParticleUpdate(SDat);
	    
	      Sim->ptrScheduler->fullUpdate(p1, p2);
	    
	      BOOST_FOREACH(shared_ptr<OutputPlugin>& Ptr, 
			     Sim->outputPlugins)
		Ptr->eventUpdate(*this, SDat, 0.0);
	    }
	}
    }
  }

  void
  SysRingDSMC::initialise(size_t nID)
  {
    ID = nID;
    dt = tstep;

    n12 = 0;
    n13 = 0;

    factor12 = range1->size()
      * diameter * M_PI * chi12 * tstep 
      / Sim->getSimVolume();
  
    factor13 = range1->size()
      * diameter * M_PI * chi13 * tstep 
      / Sim->getSimVolume();
  
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
	    Particle& p1(Sim->particles[*(range1->begin() + 2 * pairID)]);
	    Particle& p2(Sim->particles[*(range1->begin() + 2 * pairID + 1)]);
	  
	    Sim->dynamics->updateParticlePair(p1, p2);
	  
	    Vector rij;
	    for (size_t iDim(0); iDim < NDIM; ++iDim)
	      rij[iDim] = Sim->normal_sampler();
	    rij *= diameter / rij.nrm();
	    Sim->dynamics->DSMCSpheresTest(p1, p2, maxprob12, factor12, rij);
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
	    Particle& p1(Sim->particles[*(range1->begin() + id1sampler())]);

	    size_t secondID(id1sampler());

	    while ((secondID == p1.getID())
		   || ((secondID % 2) 
		       ? ((secondID-1) == p1.getID())
		       : ((secondID+1) == p1.getID())))
	      secondID = id1sampler();

	    Particle& p2(Sim->particles[*(range1->begin() + secondID)]);
	  
	    Sim->dynamics->updateParticlePair(p1, p2);
	  
	    Vector rij;
	    for (size_t iDim(0); iDim < NDIM; ++iDim)
	      rij[iDim] = Sim->normal_sampler();
	    rij *= diameter / rij.nrm();
	  
	    Sim->dynamics->DSMCSpheresTest(p1, p2, maxprob13, factor13, rij);
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
  SysRingDSMC::operator<<(const magnet::xml::Node& XML)
  {
    if (strcmp(XML.getAttribute("Type"),"RingDSMC"))
      M_throw() << "Attempting to load RingDSMC from a " 
		<< XML.getAttribute("Type") <<  " entry"; 
  
    try {
      tstep = XML.getAttribute("tStep").as<double>() * Sim->units.unitTime();    
      chi12 = XML.getAttribute("Chi12").as<double>();
      chi13 = XML.getAttribute("Chi13").as<double>();
      sysName = XML.getAttribute("Name");
      diameter = XML.getAttribute("Diameter").as<double>() * Sim->units.unitLength();
      e = XML.getAttribute("Inelasticity").as<double>();
      d2 = diameter * diameter;
      range1 = shared_ptr<Range>(Range::getClass(XML.getNode("Range1"), Sim));

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
  SysRingDSMC::outputXML(magnet::xml::XmlStream& XML) const
  {
    if (n12 || n13)
      dout << "Number of T(1,2) events " << n12
	   << "\nNumber of T(1,3) events " << n13
	   << "\nRatio T(1,2)/total " << ((double) n12) / (((double) n13) + ((double) n12))
	   << std::endl;

    XML << magnet::xml::tag("System")
	<< magnet::xml::attr("Type") << "RingDSMC"
	<< magnet::xml::attr("tStep") << tstep / Sim->units.unitTime()
	<< magnet::xml::attr("Chi12") << chi12
	<< magnet::xml::attr("Chi13") << chi13
	<< magnet::xml::attr("Diameter") << diameter / Sim->units.unitLength()
	<< magnet::xml::attr("Inelasticity") << e
	<< magnet::xml::attr("Name") << sysName
	<< magnet::xml::attr("MaxProbability12") << maxprob12
	<< magnet::xml::attr("MaxProbability13") << maxprob13
	<< magnet::xml::tag("Range1")
	<< range1
	<< magnet::xml::endtag("Range1")
	<< magnet::xml::endtag("System");
  }
}
