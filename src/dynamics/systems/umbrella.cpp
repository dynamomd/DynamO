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

#include "umbrella.hpp"
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

CSUmbrella::CSUmbrella(const XMLNode& XML, DYNAMO::SimData* tmp): 
  CSystem(tmp),
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

CSUmbrella::CSUmbrella(DYNAMO::SimData* nSim, Iflt na, Iflt nb, Iflt ndelu, 
		       std::string nName, CRange* r1, CRange* r2):
  CSystem(nSim),
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
  Iflt locdt = dt;
  
#ifdef DYNAMO_DEBUG 
  if (isnan(locdt))
    D_throw() << "A NAN system event time has been found";
#endif
  
  Sim->dSysTime += locdt;
  
  Sim->ptrScheduler->stream(locdt);
  
  //dynamics must be updated first
  Sim->Dynamics.stream(locdt);

  ++Sim->lNColl;

  BOOST_FOREACH(const size_t& id, *range1)
    Sim->Dynamics.Liouvillean().updateParticle(Sim->vParticleList[id]);
  
  BOOST_FOREACH(const size_t& id, *range2)
    Sim->Dynamics.Liouvillean().updateParticle(Sim->vParticleList[id]);
  
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

  CNParticleData SDat(Sim->Dynamics.Liouvillean().multibdyWellEvent
		      (*range1, *range2, 0.0, (kedown) ? -delU : delU, etype));

  I_cerr() << "ulevel was " << ulevel;
  if (etype != BOUNCE)
    ulevel = newulevel;

  I_cerr() << "ulevel is " << ulevel;

  Sim->signalParticleUpdate(SDat);
  
  //Only 1ParticleEvents occur
  BOOST_FOREACH(const C1ParticleData& PDat, SDat.L1partChanges)
    Sim->ptrScheduler->fullUpdate(PDat.getParticle());
  
  BOOST_FOREACH(smrtPlugPtr<COutputPlugin>& Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(*this, CNParticleData(), locdt); 
}

void
CSUmbrella::initialise(size_t nID)
{
  ID = nID;

  BOOST_FOREACH(const size_t& id, *range1)
    Sim->Dynamics.Liouvillean().updateParticle(Sim->vParticleList[id]);
  
  BOOST_FOREACH(const size_t& id, *range2)
    Sim->Dynamics.Liouvillean().updateParticle(Sim->vParticleList[id]);
  
  CPDData partdata(*Sim, *range1, *range2);

  ulevelcenter = int( - a * b * b / delU);

  Iflt r = partdata.rij.length();

  if (!ulevelset)
    {
      ulevel = int(a * (r - b) * (r - b) / delU);
      if (r < b) ulevel *= -1;
      ulevelset = true;
    }
  
  recalculateTime();

  Sim->registerParticleUpdateFunc
    (fastdelegate::MakeDelegate(this, &CSUmbrella::particlesUpdated));
}

void 
CSUmbrella::recalculateTime()
{
  BOOST_FOREACH(const size_t& id, *range1)
    Sim->Dynamics.Liouvillean().updateParticle(Sim->vParticleList[id]);
  
  BOOST_FOREACH(const size_t& id, *range2)
    Sim->Dynamics.Liouvillean().updateParticle(Sim->vParticleList[id]);
  
  CPDData partdata(*Sim, *range1, *range2);

  Iflt R_max, R_min;

  dt = HUGE_VAL;
  type = NONE;

  if (ulevel == ulevelcenter)
    {
      R_max = b - sqrt((ulevel * delU) / a);      
    
      if (b==0)//Allow a double width well if b==0
	R_max = b + sqrt((ulevel + 1 * delU) / a);
      
      //Just look for escaping as we're in the well step spanning r = 0 
      if (Sim->Dynamics.Liouvillean().SphereSphereOutRoot
	  (partdata, R_max * R_max))
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

  if (Sim->Dynamics.Liouvillean().SphereSphereInRoot(partdata, R_min * R_min))
    {
      dt = partdata.dt;
      type = WELL_IN;
    }
  else 
    if (Sim->Dynamics.Liouvillean().SphereSphereOutRoot(partdata, R_max * R_max))
      {
	dt = partdata.dt;
	type = WELL_OUT;
      }
}

void 
CSUmbrella::particlesUpdated(const CNParticleData& PDat)
{
  BOOST_FOREACH(const C1ParticleData& pdat, PDat.L1partChanges)
    if (range1->isInRange(pdat.getParticle())
	|| range2->isInRange(pdat.getParticle()))
      {
	recalculateTime();
	Sim->ptrScheduler->rebuildSystemEvents();
      }

  BOOST_FOREACH(const C2ParticleData& pdat, PDat.L2partChanges)
    if (range1->isInRange(pdat.particle1_.getParticle())
	|| range2->isInRange(pdat.particle1_.getParticle())
	|| range1->isInRange(pdat.particle2_.getParticle())
	|| range2->isInRange(pdat.particle2_.getParticle()))
      {
	recalculateTime();
	Sim->ptrScheduler->rebuildSystemEvents();
      }
}

void
CSUmbrella::operator<<(const XMLNode& XML)
{
  if (strcmp(XML.getAttribute("Type"),"Umbrella"))
    D_throw() << "Attempting to load Umbrella from a " 
	      << XML.getAttribute("Type") <<  " entry"; 
  
  try {
    sysName = XML.getAttribute("Name");

    a = boost::lexical_cast<Iflt>(XML.getAttribute("a"))
      * Sim->Dynamics.units().unitEnergy() 
      / Sim->Dynamics.units().unitArea();

    b = boost::lexical_cast<Iflt>(XML.getAttribute("b"))
      * Sim->Dynamics.units().unitLength();

    delU = boost::lexical_cast<Iflt>(XML.getAttribute("delU"))
      * Sim->Dynamics.units().unitEnergy();

    range1.set_ptr(CRange::loadClass(XML.getChildNode("Range1"), Sim));

    range2.set_ptr(CRange::loadClass(XML.getChildNode("Range2"), Sim));
    
    if (XML.isAttributeSet("currentulevel"))
      {
	ulevel = boost::lexical_cast<size_t>(XML.getAttribute("currentulevel"));
	ulevelset = true;
      }
    
  }
  catch (boost::bad_lexical_cast &)
    { D_throw() << "Failed a lexical cast in CSUmbrella"; }
}

void 
CSUmbrella::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::tag("System")
      << xmlw::attr("Type") << "Umbrella"
      << xmlw::attr("a") << a * Sim->Dynamics.units().unitArea() 
    / Sim->Dynamics.units().unitEnergy()
      << xmlw::attr("b") << b / Sim->Dynamics.units().unitLength()
      << xmlw::attr("delU") << delU / Sim->Dynamics.units().unitEnergy()
      << xmlw::attr("currentulevel") << ulevel
      << xmlw::attr("Name") << sysName
      << xmlw::tag("Range1")
      << range1
      << xmlw::endtag("Range1")
      << xmlw::tag("Range2")
      << range2
      << xmlw::endtag("Range2")
      << xmlw::endtag("System");
}
