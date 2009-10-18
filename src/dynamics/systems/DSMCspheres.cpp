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

#include "DSMCspheres.hpp"
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random/uniform_int.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../../base/is_exception.hpp"
#include "../dynamics.hpp"
#include "../units/units.hpp"
#include "../BC/BC.hpp"
#include "../../simulation/particle.hpp"
#include "../species/species.hpp"
#include "../NparticleEventData.hpp"
#include "../ranges/include.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../../schedulers/scheduler.hpp"

CSDSMCSpheres::CSDSMCSpheres(const XMLNode& XML, DYNAMO::SimData* tmp): 
  CSystem(tmp),
  uniformRand(Sim->ranGenerator, boost::uniform_real<>(0,1)),
  maxprob(0.0),
  range1(NULL),
  range2(NULL)
{
  dt = HUGE_VAL;
  operator<<(XML);
  type = DSMC;
}

CSDSMCSpheres::CSDSMCSpheres(DYNAMO::SimData* nSim, Iflt nd, Iflt ntstp, Iflt nChi, 
			     Iflt ne, std::string nName, CRange* r1, CRange* r2):
  CSystem(nSim),
  uniformRand(Sim->ranGenerator,boost::uniform_real<>(0,1)),
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
  Iflt locdt = dt;
  
#ifdef DYNAMO_DEBUG 
  if (isnan(locdt))
    D_throw() << "A NAN system event time has been found";
#endif
    
  Sim->dSysTime += locdt;
    
  Sim->ptrScheduler->stream(locdt);
  
  //dynamics must be updated first
  Sim->dynamics.stream(locdt);

  dt = tstep;

  locdt +=  Sim->freestreamAcc;
  Sim->freestreamAcc = 0;

  boost::variate_generator
    <DYNAMO::baseRNG&, boost::uniform_int<size_t> >
    id1sampler(Sim->ranGenerator, 
	      boost::uniform_int<size_t>(0, range1->size() - 1));

  boost::variate_generator
    <DYNAMO::baseRNG&, boost::uniform_int<size_t> >
    id2sampler(Sim->ranGenerator, 
	       boost::uniform_int<size_t>(0, range2->size() - 1));

  Iflt intPart;
  Iflt fracpart = std::modf(0.5 * maxprob * range1->size(),
			    &intPart);
 
  size_t nmax = static_cast<size_t>(intPart);
  
  BOOST_FOREACH(smrtPlugPtr<COutputPlugin>& Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(*this, CNParticleData(), locdt);

  if (Sim->uniform_sampler() < fracpart)
    ++nmax;

  for (size_t n = 0; n < nmax; ++n)
    {
      const CParticle& p1(Sim->vParticleList[*(range1->begin() + id1sampler())]);
      
      size_t p2id = *(range2->begin() + id2sampler());
      
      while (p2id == p1.getID())
	p2id = *(range2->begin()+id2sampler());
      
      const CParticle& p2(Sim->vParticleList[p2id]);
      
      Sim->dynamics.getLiouvillean().updateParticlePair(p1, p2);
      
      CPDData PDat;
      
      for (size_t iDim(0); iDim < NDIM; ++iDim)
	PDat.rij[iDim] = Sim->normal_sampler();
      
      PDat.rij *= diameter / PDat.rij.nrm();
      
      if (Sim->dynamics.getLiouvillean().DSMCSpheresTest
	  (p1, p2, maxprob, factor, PDat))
	{
	  ++Sim->lNColl;
	 
	  const C2ParticleData
	    SDat(Sim->dynamics.getLiouvillean().DSMCSpheresRun(p1, p2, e, PDat));

	  Sim->signalParticleUpdate(SDat);
  
	  Sim->ptrScheduler->fullUpdate(p1, p2);
	  
	  BOOST_FOREACH(smrtPlugPtr<COutputPlugin>& Ptr, Sim->outputPlugins)
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
    * diameter * PI * chi * tstep 
    / Sim->dynamics.units().simVolume();
  
  if (maxprob == 0.0)
    {
      boost::variate_generator
	<DYNAMO::baseRNG&, boost::uniform_int<size_t> >
	id1sampler(Sim->ranGenerator, 
		   boost::uniform_int<size_t>(0, range1->size() - 1));
      
      boost::variate_generator
	<DYNAMO::baseRNG&, boost::uniform_int<size_t> >
	id2sampler(Sim->ranGenerator, 
		   boost::uniform_int<size_t>(0, range2->size() - 1));
      
      //Just do some quick testing to get an estimate
      for (size_t n = 0; n < 1000; ++n)
	{
	  const CParticle& p1(Sim->vParticleList[*(range1->begin() + id1sampler())]);
	  
	  size_t p2id = *(range2->begin() + id2sampler());
	  
	  while (p2id == p1.getID())
	    p2id = *(range2->begin()+id2sampler());
	  
	  const CParticle& p2(Sim->vParticleList[p2id]);
	  
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
    I_cerr() << "MaxProbability is " << maxprob
	     << "\nNpairs per step is " << 0.5 * range1->size() * maxprob;
  else
    I_cout() << "MaxProbability is " << maxprob
	     << "\nNpairs per step is " << 0.5 * range1->size() * maxprob;
  
  if (0.5 * range1->size() * maxprob < 2.0)
    I_cerr() << "This probability is low";
}

void
CSDSMCSpheres::operator<<(const XMLNode& XML)
{
  if (strcmp(XML.getAttribute("Type"),"DSMCSpheres"))
    D_throw() << "Attempting to load DSMCSpheres from a " << XML.getAttribute("Type") <<  " entry"; 
  
  try {
    tstep = boost::lexical_cast<Iflt>(XML.getAttribute("tStep"))
      * Sim->dynamics.units().unitTime();
    
    chi = boost::lexical_cast<Iflt>(XML.getAttribute("Chi"));
    
    sysName = XML.getAttribute("Name");

    diameter = boost::lexical_cast<Iflt>(XML.getAttribute("Diameter"))
      * Sim->dynamics.units().unitLength();

    e = boost::lexical_cast<Iflt>(XML.getAttribute("Inelasticity"));

    d2 = diameter * diameter;

    range1.set_ptr(CRange::loadClass(XML.getChildNode("Range1"), Sim));

    range2.set_ptr(CRange::loadClass(XML.getChildNode("Range2"), Sim));

    if (XML.isAttributeSet("MaxProbability"))
      maxprob = boost::lexical_cast<Iflt>(XML.getAttribute("MaxProbability"));
	
  }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in CGGlobal";
    }
}

void 
CSDSMCSpheres::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::tag("System")
      << xmlw::attr("Type") << "DSMCSpheres"
      << xmlw::attr("tStep") << tstep / Sim->dynamics.units().unitTime()
      << xmlw::attr("Chi") << chi
      << xmlw::attr("Diameter") << diameter / Sim->dynamics.units().unitLength()
      << xmlw::attr("Inelasticity") << e
      << xmlw::attr("Name") << sysName
      << xmlw::attr("MaxProbability") << maxprob
      << xmlw::tag("Range1")
      << range1
      << xmlw::endtag("Range1")
      << xmlw::tag("Range2")
      << range2
      << xmlw::endtag("Range2")
      << xmlw::endtag("System");
}
