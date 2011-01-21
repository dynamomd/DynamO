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

#include "umbrella.hpp"
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random/uniform_int.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../dynamics.hpp"
#include "../units/units.hpp"
#include "../BC/BC.hpp"
#include "../../simulation/particle.hpp"
#include "../species/species.hpp"
#include "../NparticleEventData.hpp"
#include "../ranges/include.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../../schedulers/scheduler.hpp"

CSUmbrella::CSUmbrella(const XMLNode& XML, DYNAMO::SimData* tmp): 
  System(tmp),
  a(1.0),
  b(1.0),
  delU(0.1),
  ulevelcenter(0),
  ulevel(-1),
  ulevelset(false),
  range1(NULL),
  range2(NULL)
{
  dt = HUGE_VAL;
  operator<<(XML);
  type = UMBRELLA;
}

CSUmbrella::CSUmbrella(DYNAMO::SimData* nSim, double na, double nb, double ndelu, 
		       std::string nName, CRange* r1, CRange* r2):
  System(nSim),
  a(na),
  b(nb),
  delU(ndelu),
  ulevelcenter(0),
  ulevel(-1),
  ulevelset(false),
  range1(r1),
  range2(r2)
{
  sysName = nName;
  type = UMBRELLA;
}

void 
CSUmbrella::runEvent() const
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

  ++Sim->eventCount;

  BOOST_FOREACH(const size_t& id, *range1)
    Sim->dynamics.getLiouvillean().updateParticle(Sim->particleList[id]);
  
  BOOST_FOREACH(const size_t& id, *range2)
    Sim->dynamics.getLiouvillean().updateParticle(Sim->particleList[id]);
  
  bool kedown(false); //Will kinetic energy go down?

  int newulevel ;

  if (ulevel == 0)
    {
      kedown = true;
      
      if (type == WELL_OUT)
	newulevel = 1;
      else
	newulevel = -1;
    }
  else if (type == WELL_OUT)
    {
      if (ulevel > 0) kedown = true;
      newulevel = ulevel + 1; 
    }
  else //if (type == WELL_IN)
    {
      if (ulevel < 0) kedown = true;
      newulevel = ulevel - 1;
    }
    
  EEventType etype(NONE);

  NEventData SDat(Sim->dynamics.getLiouvillean().multibdyWellEvent
		      (*range1, *range2, 0.0, (kedown) ? -delU : delU, etype));

  if (etype != BOUNCE)
    ulevel = newulevel;

  Sim->signalParticleUpdate(SDat);
  
  //Only 1ParticleEvents occur
  BOOST_FOREACH(const ParticleEventData& PDat, SDat.L1partChanges)
    Sim->ptrScheduler->fullUpdate(PDat.getParticle());
  
  locdt += Sim->freestreamAcc;

  Sim->freestreamAcc = 0;

  BOOST_FOREACH(magnet::ClonePtr<OutputPlugin>& Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(*this, SDat, locdt); 
}

void
CSUmbrella::initialise(size_t nID)
{
  ID = nID;

  BOOST_FOREACH(const size_t& id, *range1)
    Sim->dynamics.getLiouvillean().updateParticle(Sim->particleList[id]);
  
  BOOST_FOREACH(const size_t& id, *range2)
    Sim->dynamics.getLiouvillean().updateParticle(Sim->particleList[id]);
  
  CPDData partdata(*Sim, *range1, *range2);

  ulevelcenter = int( - a * b * b / delU);

  double r = partdata.rij.nrm();

  if (!ulevelset)
    {
      ulevel = int(a * (r - b) * (r - b) / delU);
      if (r < b) ulevel *= -1;
      ulevelset = true;
    }
  
  recalculateTime();

  Sim->registerParticleUpdateFunc
    (magnet::function::MakeDelegate(this, &CSUmbrella::particlesUpdated));
}

void 
CSUmbrella::recalculateTime()
{
  BOOST_FOREACH(const size_t& id, *range1)
    Sim->dynamics.getLiouvillean().updateParticle(Sim->particleList[id]);
  
  BOOST_FOREACH(const size_t& id, *range2)
    Sim->dynamics.getLiouvillean().updateParticle(Sim->particleList[id]);
  
  CPDData partdata(*Sim, *range1, *range2);

  double R_max, R_min;

  dt = HUGE_VAL;
  type = NONE;

  if (ulevel == ulevelcenter)
    {
      R_max = b - sqrt((ulevel * delU) / a);      
    
      if (b==0)//Allow a double width well if b==0
	R_max = b + sqrt((ulevel + 1 * delU) / a);
      
      //Just look for escaping as we're in the well step spanning r = 0 
      if (Sim->dynamics.getLiouvillean().SphereSphereOutRoot
	  (partdata, R_max * R_max, true, true))
	{
	  dt = partdata.dt;
	  type = WELL_OUT;
	}
      
      return;
    }
  
  if (ulevel == 0)
    {
      //We're on the minimum

      //We don't worry about the minimum crossing r=0, as this is
      //caught by the above if statement
      
      R_max = b + sqrt((1 * delU) / a);
      R_min = b - sqrt((1 * delU) / a);
    }
  else if (ulevel < 0)
    {
      R_max = b - sqrt((-ulevel) * delU / a);
      R_min = b - sqrt(((-ulevel) + 1) * delU / a);
    }
  else
    {
      R_min = b + sqrt((ulevel * delU) / a);
      R_max = b + sqrt(((ulevel + 1) * delU) / a);
    }

  if (Sim->dynamics.getLiouvillean().SphereSphereInRoot(partdata, R_min * R_min, true, true))
    {
      dt = partdata.dt;
      type = WELL_IN;
    }
  else 
    if (Sim->dynamics.getLiouvillean().SphereSphereOutRoot(partdata, R_max * R_max, true, true))
      {
	dt = partdata.dt;
	type = WELL_OUT;
      }
}

void 
CSUmbrella::particlesUpdated(const NEventData& PDat)
{
  BOOST_FOREACH(const ParticleEventData& pdat, PDat.L1partChanges)
    if (range1->isInRange(pdat.getParticle())
	|| range2->isInRange(pdat.getParticle()))
      {
	recalculateTime();
	Sim->ptrScheduler->rebuildSystemEvents();
	return;
      }

  BOOST_FOREACH(const PairEventData& pdat, PDat.L2partChanges)
    if (range1->isInRange(pdat.particle1_.getParticle())
	|| range2->isInRange(pdat.particle1_.getParticle())
	|| range1->isInRange(pdat.particle2_.getParticle())
	|| range2->isInRange(pdat.particle2_.getParticle()))
      {
	recalculateTime();
	Sim->ptrScheduler->rebuildSystemEvents();
	return;
      }
}

void
CSUmbrella::operator<<(const XMLNode& XML)
{
  if (strcmp(XML.getAttribute("Type"),"Umbrella"))
    M_throw() << "Attempting to load Umbrella from a " 
	      << XML.getAttribute("Type") <<  " entry"; 
  
  try {
    sysName = XML.getAttribute("Name");

    a = boost::lexical_cast<double>(XML.getAttribute("a"))
      * Sim->dynamics.units().unitEnergy() 
      / Sim->dynamics.units().unitArea();

    b = boost::lexical_cast<double>(XML.getAttribute("b"))
      * Sim->dynamics.units().unitLength();

    delU = boost::lexical_cast<double>(XML.getAttribute("delU"))
      * Sim->dynamics.units().unitEnergy();

    range1.set_ptr(CRange::loadClass(XML.getChildNode("Range1"), Sim));

    range2.set_ptr(CRange::loadClass(XML.getChildNode("Range2"), Sim));
    
    if (XML.isAttributeSet("currentulevel"))
      {
	ulevel = boost::lexical_cast<size_t>(XML.getAttribute("currentulevel"));
	ulevelset = true;
      }
    
  }
  catch (boost::bad_lexical_cast &)
    { M_throw() << "Failed a lexical cast in CSUmbrella"; }
}

void 
CSUmbrella::outputXML(xml::XmlStream& XML) const
{
  XML << xml::tag("System")
      << xml::attr("Type") << "Umbrella"
      << xml::attr("a") << a * Sim->dynamics.units().unitArea() 
    / Sim->dynamics.units().unitEnergy()
      << xml::attr("b") << b / Sim->dynamics.units().unitLength()
      << xml::attr("delU") << delU / Sim->dynamics.units().unitEnergy()
      << xml::attr("currentulevel") << ulevel
      << xml::attr("Name") << sysName
      << xml::tag("Range1")
      << range1
      << xml::endtag("Range1")
      << xml::tag("Range2")
      << range2
      << xml::endtag("Range2")
      << xml::endtag("System");
}
