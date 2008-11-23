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
  range(NULL)
{
  dt = HUGE_VAL;
  operator<<(XML);
  type = DSMC;
}

CSDSMCSpheres::CSDSMCSpheres(DYNAMO::SimData* nSim, Iflt nd, Iflt ntstp, Iflt nChi, 
			     Iflt ne,
			     std::string nName):
  CSystem(nSim),
  uniformRand(Sim->ranGenerator,boost::uniform_real<>(0,1)),
  tstep(ntstp),
  chi(nChi),
  d2(nd * nd),
  diameter(nd),
  maxprob(0.0),
  e(ne),
  range(new CRAll(Sim))
{
  sysName = nName;
  type = DSMC;
}

void 
CSDSMCSpheres::stream(Iflt ndt)
{
  dt -= ndt;
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
  Sim->Dynamics.stream(locdt);

  dt = tstep;

  CNParticleData SDat;

  boost::variate_generator
    <DYNAMO::baseRNG&, boost::uniform_int<unsigned int> >
    idsampler(Sim->ranGenerator, 
	      boost::uniform_int<unsigned int>(0, range->size() - 1));

  Iflt intPart;
  Iflt fracpart = std::modf(0.5 * maxprob * range->size(), 
			    &intPart);
 
  size_t nmax = static_cast<size_t>(intPart);

  if (Sim->uniform_sampler() < fracpart)
    ++nmax;

  for (size_t n = 0; n < nmax; ++n)
    {
      const CParticle& p1(Sim->vParticleList[*(range->begin() + idsampler())]);
      
      size_t p2id = *(range->begin() + idsampler());
      
      while (p2id == p1.getID())
	p2id = *(range->begin()+idsampler());
      
      const CParticle& p2(Sim->vParticleList[p2id]);
      
      Sim->Dynamics.Liouvillean().updateParticlePair(p1, p2);
      
      CPDData PDat;
      
      for (size_t iDim(0); iDim < NDIM; ++iDim)
	PDat.rij[iDim] = Sim->normal_sampler();
      
      PDat.rij *= diameter / PDat.rij.length();
      
      if (Sim->Dynamics.Liouvillean().DSMCSpheresTest
	  (p1, p2, maxprob, factor, PDat))
	{
	  ++Sim->lNColl;
	  
	  SDat.L2partChanges.push_back
	    (Sim->Dynamics.Liouvillean().DSMCSpheresRun(p1, p2, e, PDat));
	}
    }
    
  BOOST_FOREACH(smrtPlugPtr<COutputPlugin>& Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(*this, SDat, locdt);
}

void
CSDSMCSpheres::initialise(size_t nID)
{
  ID = nID;
  dt = tstep;

  factor = 4.0 * range->size()
    * diameter * PI * chi * tstep 
    / Sim->Dynamics.units().simVolume();
  
  if (maxprob == 0.0)
    BOOST_FOREACH(const size_t& id1, *range)
      BOOST_FOREACH(const size_t& id2, *range)
      // Test everything to get an estimate for the max probibility
      if (id1 != id2)
	{
	  const CParticle& p1(Sim->vParticleList[id1]);
	  const CParticle& p2(Sim->vParticleList[id2]);

	  Sim->Dynamics.Liouvillean().updateParticlePair(p1, p2);

	  CPDData PDat;
	  
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    PDat.rij[iDim] = Sim->normal_sampler();
	  
	  PDat.rij *= diameter / PDat.rij.length();
	  
	  Sim->Dynamics.Liouvillean().DSMCSpheresTest(p1, p2, maxprob, 
						      factor, PDat);
	}
}

void
CSDSMCSpheres::operator<<(const XMLNode& XML)
{
  if (strcmp(XML.getAttribute("Type"),"DSMCSpheres"))
    D_throw() << "Attempting to load DSMCSpheres from a " << XML.getAttribute("Type") <<  " entry"; 
  
  try {
    tstep = boost::lexical_cast<Iflt>(XML.getAttribute("tStep"))
      * Sim->Dynamics.units().unitTime();
    
    chi = boost::lexical_cast<Iflt>(XML.getAttribute("Chi"));
    
    sysName = XML.getAttribute("Name");

    diameter = boost::lexical_cast<Iflt>(XML.getAttribute("Diameter"))
      * Sim->Dynamics.units().unitLength();

    e = boost::lexical_cast<Iflt>(XML.getAttribute("Inelasticity"));

    d2 = diameter * diameter;

    range.set_ptr(CRange::loadClass(XML,Sim));

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
      << xmlw::attr("tStep") << tstep / Sim->Dynamics.units().unitTime()
      << xmlw::attr("Chi") << chi
      << xmlw::attr("Diameter") << diameter / Sim->Dynamics.units().unitLength()
      << xmlw::attr("Inelasticity") << e
      << xmlw::attr("Name") << sysName
      << xmlw::attr("MaxProbability") << maxprob
      << range
      << xmlw::endtag("System");
}
