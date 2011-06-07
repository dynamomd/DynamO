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

#include "DSMCspheres.hpp"
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

CSDSMCSpheres::CSDSMCSpheres(const magnet::xml::Node& XML, dynamo::SimData* tmp): 
  System(tmp),
  maxprob(0.0),
  range1(NULL),
  range2(NULL)
{
  dt = HUGE_VAL;
  operator<<(XML);
  type = DSMC;
}

CSDSMCSpheres::CSDSMCSpheres(dynamo::SimData* nSim, double nd, double ntstp, double nChi, 
			     double ne, std::string nName, CRange* r1, CRange* r2):
  System(nSim),
  tstep(ntstp),
  chi(nChi),
  d2(nd * nd),
  diameter(nd),
  maxprob(0.0),
  e(ne),
  range1(r1),
  range2(r2)
{
  sysName = nName;
  type = DSMC;
}

void 
CSDSMCSpheres::runEvent() const
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

  locdt +=  Sim->freestreamAcc;
  Sim->freestreamAcc = 0;

  boost::variate_generator
    <dynamo::baseRNG&, boost::uniform_int<size_t> >
    id1sampler(Sim->ranGenerator, 
	      boost::uniform_int<size_t>(0, range1->size() - 1));

  boost::variate_generator
    <dynamo::baseRNG&, boost::uniform_int<size_t> >
    id2sampler(Sim->ranGenerator, 
	       boost::uniform_int<size_t>(0, range2->size() - 1));

  double intPart;
  double fracpart = std::modf(0.5 * maxprob * range1->size(),
			    &intPart);
 
  size_t nmax = static_cast<size_t>(intPart);
  
  BOOST_FOREACH(magnet::ClonePtr<OutputPlugin>& Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(*this, NEventData(), locdt);

  if (Sim->uniform_sampler() < fracpart)
    ++nmax;

  for (size_t n = 0; n < nmax; ++n)
    {
      const Particle& p1(Sim->particleList[*(range1->begin() + id1sampler())]);
      
      size_t p2id = *(range2->begin() + id2sampler());
      
      while (p2id == p1.getID())
	p2id = *(range2->begin()+id2sampler());
      
      const Particle& p2(Sim->particleList[p2id]);
      
      Sim->dynamics.getLiouvillean().updateParticlePair(p1, p2);
      
      CPDData PDat;
      
      for (size_t iDim(0); iDim < NDIM; ++iDim)
	PDat.rij[iDim] = Sim->normal_sampler();
      
      PDat.rij *= diameter / PDat.rij.nrm();
      
      if (Sim->dynamics.getLiouvillean().DSMCSpheresTest
	  (p1, p2, maxprob, factor, PDat))
	{
	  ++Sim->eventCount;
	 
	  const PairEventData
	    SDat(Sim->dynamics.getLiouvillean().DSMCSpheresRun(p1, p2, e, PDat));

	  Sim->signalParticleUpdate(SDat);
  
	  Sim->ptrScheduler->fullUpdate(p1, p2);
	  
	  BOOST_FOREACH(magnet::ClonePtr<OutputPlugin>& Ptr, Sim->outputPlugins)
	    Ptr->eventUpdate(*this, SDat, 0.0);
	}
    }

}

void
CSDSMCSpheres::initialise(size_t nID)
{
  ID = nID;
  dt = tstep;

  factor = 4.0 * range2->size()
    * diameter * M_PI * chi * tstep 
    / Sim->dynamics.getSimVolume();
  
  if (maxprob == 0.0)
    {
      boost::variate_generator
	<dynamo::baseRNG&, boost::uniform_int<size_t> >
	id1sampler(Sim->ranGenerator, 
		   boost::uniform_int<size_t>(0, range1->size() - 1));
      
      boost::variate_generator
	<dynamo::baseRNG&, boost::uniform_int<size_t> >
	id2sampler(Sim->ranGenerator, 
		   boost::uniform_int<size_t>(0, range2->size() - 1));
      
      //Just do some quick testing to get an estimate
      for (size_t n = 0; n < 1000; ++n)
	{
	  const Particle& p1(Sim->particleList[*(range1->begin() + id1sampler())]);
	  
	  size_t p2id = *(range2->begin() + id2sampler());
	  
	  while (p2id == p1.getID())
	    p2id = *(range2->begin()+id2sampler());
	  
	  const Particle& p2(Sim->particleList[p2id]);
	  
	  Sim->dynamics.getLiouvillean().updateParticlePair(p1, p2);
	  
	  CPDData PDat;
	  
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    PDat.rij[iDim] = Sim->normal_sampler();
	
	  PDat.rij *= diameter / PDat.rij.nrm();
	  
	  Sim->dynamics.getLiouvillean().DSMCSpheresTest(p1, p2, maxprob, 
						      factor, PDat);
	}
    }

  if (maxprob > 0.5)
    derr << "MaxProbability is " << maxprob
	     << "\nNpairs per step is " << 0.5 * range1->size() * maxprob << std::endl;
  else
    dout << "MaxProbability is " << maxprob
	     << "\nNpairs per step is " << 0.5 * range1->size() * maxprob << std::endl;
  
  if (0.5 * range1->size() * maxprob < 2.0)
    derr << "This probability is low" << std::endl;
}

void
CSDSMCSpheres::operator<<(const magnet::xml::Node& XML)
{
  if (strcmp(XML.getAttribute("Type"),"DSMCSpheres"))
    M_throw() << "Attempting to load DSMCSpheres from a " << XML.getAttribute("Type") <<  " entry"; 
  
  try {
    tstep = XML.getAttribute("tStep").as<double>() * Sim->dynamics.units().unitTime();
    chi = XML.getAttribute("Chi").as<double>();
    sysName = XML.getAttribute("Name");
    diameter = XML.getAttribute("Diameter").as<double>() * Sim->dynamics.units().unitLength();
    e = XML.getAttribute("Inelasticity").as<double>();
    d2 = diameter * diameter;
    range1.set_ptr(CRange::getClass(XML.getNode("Range1"), Sim));
    range2.set_ptr(CRange::getClass(XML.getNode("Range2"), Sim));
    if (XML.getAttribute("MaxProbability").valid())
      maxprob = XML.getAttribute("MaxProbability").as<double>();
  }
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in CGGlobal";
    }
}

void 
CSDSMCSpheres::outputXML(xml::XmlStream& XML) const
{
  XML << xml::tag("System")
      << xml::attr("Type") << "DSMCSpheres"
      << xml::attr("tStep") << tstep / Sim->dynamics.units().unitTime()
      << xml::attr("Chi") << chi
      << xml::attr("Diameter") << diameter / Sim->dynamics.units().unitLength()
      << xml::attr("Inelasticity") << e
      << xml::attr("Name") << sysName
      << xml::attr("MaxProbability") << maxprob
      << xml::tag("Range1")
      << range1
      << xml::endtag("Range1")
      << xml::tag("Range2")
      << range2
      << xml::endtag("Range2")
      << xml::endtag("System");
}
