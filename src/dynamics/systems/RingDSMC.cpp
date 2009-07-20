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

#include "RingDSMC.hpp"
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

CSRingDSMC::CSRingDSMC(const XMLNode& XML, DYNAMO::SimData* tmp): 
  CSystem(tmp),
  uniformRand(Sim->ranGenerator, boost::uniform_real<>(0,1)),
  maxprob12(0.0),
  maxprob13(0.0),
  range1(NULL)
{
  dt = HUGE_VAL;
  operator<<(XML);
  type = DSMC;
}

CSRingDSMC::CSRingDSMC(DYNAMO::SimData* nSim, Iflt nd, Iflt ntstp, Iflt nChi, 
			     Iflt ne, std::string nName, CRange* r1):
  CSystem(nSim),
  uniformRand(Sim->ranGenerator,boost::uniform_real<>(0,1)),
  tstep(ntstp),
  chi(nChi),
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
    D_throw() << "Need an even number of particles in"
	      << " the range to make a whole number of velocity pairs";
}

void 
CSRingDSMC::runEvent() const
{
  Iflt locdt = dt;
  
#ifdef DYNAMO_DEBUG 
  if (isnan(locdt))
    D_throw() << "A NAN system event time has been found";
#endif
    
  Sim->dSysTime += locdt;
    
  Sim->ptrScheduler->stream(locdt);
  
  //dynamics must be updated first
  Sim->Dynamics.stream(locdt);

  dt = tstep;

  locdt += Sim->freestreamAcc;
  Sim->freestreamAcc = 0;

  Iflt intPart;
  Iflt fracpart = std::modf(maxprob12 * (range1->size()/2),
			    &intPart);
 
  size_t nmax = static_cast<size_t>(intPart);
  
  if (Sim->uniform_sampler() < fracpart) ++nmax;

  BOOST_FOREACH(smrtPlugPtr<COutputPlugin>& Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(*this, CNParticleData(), locdt);
  
  CNParticleData ndata;

  boost::variate_generator
    <DYNAMO::baseRNG&, boost::uniform_int<size_t> >
    id1sampler(Sim->ranGenerator, 
	       boost::uniform_int<size_t>(0, (range1->size()/2) - 1));

  for (size_t n = 0; n < nmax; ++n)
    {
      size_t pairID(id1sampler());
      const CParticle& p1(Sim->vParticleList[*(range1->begin() + 2 * pairID)]);
      const CParticle& p2(Sim->vParticleList[*(range1->begin() + 2 * pairID + 1)]);
            
      Sim->Dynamics.Liouvillean().updateParticlePair(p1, p2);
      
      CPDData PDat;
      
      for (size_t iDim(0); iDim < NDIM; ++iDim)
	PDat.rij[iDim] = Sim->normal_sampler();
      
      PDat.rij *= diameter / PDat.rij.nrm();
      
      if (Sim->Dynamics.Liouvillean().DSMCSpheresTest
	  (p1, p2, maxprob12, factor12, PDat))
	{
	  ++Sim->lNColl;
	 
	  const C2ParticleData
	    SDat(Sim->Dynamics.Liouvillean().DSMCSpheresRun(p1, p2, e, PDat));

	  Sim->signalParticleUpdate(SDat);
  
	  Sim->ptrScheduler->fullUpdate(p1, p2);

	  ndata.L2partChanges.push_back(SDat);
	  
	  BOOST_FOREACH(smrtPlugPtr<COutputPlugin>& Ptr, Sim->outputPlugins)
	    Ptr->eventUpdate(*this, SDat, 0.0);
	}
    }

}

void
CSRingDSMC::initialise(size_t nID)
{
  ID = nID;
  dt = tstep;

  factor12 = 4.0 * range1->size()
    * diameter * PI * chi * tstep 
    / Sim->Dynamics.units().simVolume();
  
  if (maxprob12 == 0.0)
    { 
      boost::variate_generator
	<DYNAMO::baseRNG&, boost::uniform_int<size_t> >
	id1sampler(Sim->ranGenerator, 
		   boost::uniform_int<size_t>(0, (range1->size()/2) - 1));
      
      //Just do some quick testing to get an estimate
      for (size_t n = 0; n < 1000; ++n)
	{
	  size_t pairID(id1sampler());
	  const CParticle& p1(Sim->vParticleList[*(range1->begin() + 2 * pairID)]);
	  const CParticle& p2(Sim->vParticleList[*(range1->begin() + 2 * pairID + 1)]);
	  
	  Sim->Dynamics.Liouvillean().updateParticlePair(p1, p2);
	  
	  CPDData PDat;
	  
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    PDat.rij[iDim] = Sim->normal_sampler();
	  
	  PDat.rij *= diameter / PDat.rij.nrm();
	  
	  Sim->Dynamics.Liouvillean().DSMCSpheresTest(p1, p2, maxprob12, 
						      factor12, PDat);
	}
    }

  if (maxprob12 > 0.5)
    I_cerr() << "MaxProbability12 is " << maxprob12
	     << "\nNpairs per step is " << range1->size() * maxprob12;
  else
    I_cout() << "MaxProbability12 is " << maxprob12
	     << "\nNpairs per step is " << range1->size() * maxprob12;
  
  if (range1->size() * maxprob12 < 2.0)
    I_cerr() << "This probability is low";
}

void
CSRingDSMC::operator<<(const XMLNode& XML)
{
  if (strcmp(XML.getAttribute("Type"),"RingDSMC"))
    D_throw() << "Attempting to load RingDSMC from a " << XML.getAttribute("Type") <<  " entry"; 
  
  try {
    tstep = boost::lexical_cast<Iflt>(XML.getAttribute("tStep"))
      * Sim->Dynamics.units().unitTime();
    
    chi = boost::lexical_cast<Iflt>(XML.getAttribute("Chi"));
    
    sysName = XML.getAttribute("Name");

    diameter = boost::lexical_cast<Iflt>(XML.getAttribute("Diameter"))
      * Sim->Dynamics.units().unitLength();

    e = boost::lexical_cast<Iflt>(XML.getAttribute("Inelasticity"));

    d2 = diameter * diameter;

    range1.set_ptr(CRange::loadClass(XML.getChildNode("Range1"), Sim));

    if (XML.isAttributeSet("MaxProbability12"))
      maxprob12 = boost::lexical_cast<Iflt>(XML.getAttribute("MaxProbability12"));
	
    if (XML.isAttributeSet("MaxProbability13"))
      maxprob13 = boost::lexical_cast<Iflt>(XML.getAttribute("MaxProbability13"));
  }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in CGGlobal";
    }
}

void 
CSRingDSMC::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::tag("System")
      << xmlw::attr("Type") << "RingDSMC"
      << xmlw::attr("tStep") << tstep / Sim->Dynamics.units().unitTime()
      << xmlw::attr("Chi") << chi
      << xmlw::attr("Diameter") << diameter / Sim->Dynamics.units().unitLength()
      << xmlw::attr("Inelasticity") << e
      << xmlw::attr("Name") << sysName
      << xmlw::attr("MaxProbability12") << maxprob12
      << xmlw::attr("MaxProbability13") << maxprob13
      << xmlw::tag("Range1")
      << range1
      << xmlw::endtag("Range1")
      << xmlw::endtag("System");
}
