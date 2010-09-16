/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "softcore.hpp"
#include <boost/lexical_cast.hpp>
#include <cmath>
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../BC/BC.hpp"
#include "../dynamics.hpp"
#include "../units/units.hpp"
#include "../globals/global.hpp"
#include "../../simulation/particle.hpp"
#include "../interactions/intEvent.hpp"
#include "../species/species.hpp"
#include "../2particleEventData.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../../base/is_simdata.hpp"
#include <iomanip>
#include "../../schedulers/scheduler.hpp"
#include "../NparticleEventData.hpp"


ISoftCore::ISoftCore(DYNAMO::SimData* tmp, Iflt nd, Iflt nWD, 
		       C2Range* nR):
  ISingleCapture(tmp,nR),
  diameter(nd),d2(nd*nd),wellDepth(nWD) 
{}

ISoftCore::ISoftCore(const XMLNode& XML, DYNAMO::SimData* tmp):
  ISingleCapture(tmp,NULL) //A temporary value!
{
  operator<<(XML);
}

void 
ISoftCore::operator<<(const XMLNode& XML)
{
  if (strcmp(XML.getAttribute("Type"),"SoftCore"))
    M_throw() << "Attempting to load SoftCore from non SoftCore entry";
  
  range.set_ptr(C2Range::loadClass(XML,Sim));
  
  try {
    diameter = Sim->dynamics.units().unitLength() 
      * boost::lexical_cast<Iflt>(XML.getAttribute("Diameter"));
    
    wellDepth = boost::lexical_cast<Iflt>(XML.getAttribute("WellDepth"))
      * Sim->dynamics.units().unitEnergy();
    
    d2 = diameter * diameter;
    
    intName = XML.getAttribute("Name");

    ISingleCapture::loadCaptureMap(XML);   
  }
  catch (boost::bad_lexical_cast &)
    { M_throw() << "Failed a lexical cast in CISoftCore"; }
}

Interaction* 
ISoftCore::Clone() const 
{ return new ISoftCore(*this); }

Iflt 
ISoftCore::hardCoreDiam() const 
{ return diameter; }

Iflt 
ISoftCore::maxIntDist() const 
{ return diameter; }

void 
ISoftCore::rescaleLengths(Iflt scale) 
{ 
  diameter += scale*diameter; 
  d2 = diameter*diameter;
}

void 
ISoftCore::initialise(size_t nID)
{
  ID = nID;
  ISingleCapture::initCaptureMap();
}

bool 
ISoftCore::captureTest(const Particle& p1, const Particle& p2) const
{
  Vector  rij = p1.getPosition() - p2.getPosition();
  Sim->dynamics.BCs().applyBC(rij);
  
  return (rij.nrm2() <= d2);
}

IntEvent 
ISoftCore::getEvent(const Particle &p1, 
		     const Particle &p2) const 
{
#ifdef DYNAMO_DEBUG
  if (!Sim->dynamics.getLiouvillean().isUpToDate(p1))
    M_throw() << "Particle 1 is not up to date";
  
  if (!Sim->dynamics.getLiouvillean().isUpToDate(p2))
    M_throw() << "Particle 2 is not up to date";

  if (p1 == p2)
    M_throw() << "You shouldn't pass p1==p2 events to the interactions!";
#endif 

  CPDData colldat(*Sim, p1, p2);
    
  if (isCaptured(p1, p2)) 
    {
      if (Sim->dynamics.getLiouvillean().SphereSphereOutRoot(colldat, d2))
	return IntEvent(p1, p2, colldat.dt, WELL_OUT, *this);
    }
  else if (Sim->dynamics.getLiouvillean().SphereSphereInRoot(colldat, d2)) 
    {
#ifdef DYNAMO_OverlapTesting
      if (Sim->dynamics.getLiouvillean().sphereOverlap(colldat,d2))
	M_throw() << "Overlapping cores (but not registered as captured) particles found in soft core" 
		  << "\nparticle1 " << p1.getID() << ", particle2 " 
		  << p2.getID() << "\nOverlap = " 
		  << (sqrt(colldat.r2) - sqrt(d2)) / Sim->dynamics.units().unitLength();
#endif

      return IntEvent(p1, p2, colldat.dt, WELL_IN, *this);
    }

  return IntEvent(p1, p2, HUGE_VAL, NONE, *this);
}

void
ISoftCore::runEvent(const Particle& p1, const Particle& p2, const IntEvent& iEvent) const
{
  ++Sim->eventCount;
  
  switch (iEvent.getType())
    {
    case WELL_IN:
      {
	PairEventData retVal(Sim->dynamics.getLiouvillean()
			      .SphereWellEvent(iEvent, wellDepth, d2));
	
	if (retVal.getType() != BOUNCE)
	  addToCaptureMap(p1, p2);      
	
	//Now we're past the event, update the scheduler and plugins
	Sim->signalParticleUpdate(retVal);
	Sim->ptrScheduler->fullUpdate(p1, p2);
	
	BOOST_FOREACH(ClonePtr<OutputPlugin> & Ptr, Sim->outputPlugins)
	  Ptr->eventUpdate(iEvent, retVal);


	break;
      }
    case WELL_OUT:
      {
	PairEventData retVal(Sim->dynamics.getLiouvillean()
			      .SphereWellEvent(iEvent, -wellDepth, d2));
	
	if (retVal.getType() != BOUNCE)
	  removeFromCaptureMap(p1, p2);      
	
	Sim->signalParticleUpdate(retVal);

	//Now we're past the event, update the scheduler and plugins
	Sim->ptrScheduler->fullUpdate(p1, p2);
	
	BOOST_FOREACH(ClonePtr<OutputPlugin>& Ptr, 
		      Sim->outputPlugins)
	  Ptr->eventUpdate(iEvent, retVal);

	break;
      }
    default:
      M_throw() << "Unknown collision type";
    }
}

void
ISoftCore::checkOverlaps(const Particle& part1, const Particle& part2) const
{
  Vector  rij = part1.getPosition() - part2.getPosition();
  Sim->dynamics.BCs().applyBC(rij);
  Iflt r2 = rij.nrm2();

  if (isCaptured(part1, part2))
    {
      if (r2 > d2)
	I_cerr() << std::setprecision(std::numeric_limits<float>::digits10)
		 << "Possible escaped captured pair in diagnostics\n ID1=" << part1.getID() 
		 << ", ID2=" << part2.getID() << "\nR_ij^2=" 
		 << r2 / pow(Sim->dynamics.units().unitLength(),2)
		 << "\nd^2=" 
		 << d2 / pow(Sim->dynamics.units().unitLength(),2);
    }
  else 
    if (r2 < d2)
      I_cerr() << std::setprecision(std::numeric_limits<float>::digits10)
	       << "Possible missed captured pair in diagnostics\n ID1=" << part1.getID() 
	       << ", ID2=" << part2.getID() << "\nR_ij^2=" 
	       << r2 / pow(Sim->dynamics.units().unitLength(),2)
	       << "\nd^2=" 
	       << d2 / pow(Sim->dynamics.units().unitLength(),2);
}
  
void 
ISoftCore::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Type") << "SoftCore"
      << xml::attr("Diameter") 
      << diameter / Sim->dynamics.units().unitLength() 
      << xml::attr("WellDepth") 
      << wellDepth / Sim->dynamics.units().unitEnergy()
      << xml::attr("Name") << intName
      << range;
  
  ISingleCapture::outputCaptureMap(XML);  
}
