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

#include "softcore.hpp"
#include <boost/lexical_cast.hpp>
#include <cmath>
#include "../../base/is_exception.hpp"
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


CISoftCore::CISoftCore(DYNAMO::SimData* tmp, Iflt nd, Iflt nWD, 
		       C2Range* nR):
  CICapture(tmp,nR),
  diameter(nd),d2(nd*nd),wellDepth(nWD) 
{}

CISoftCore::CISoftCore(const XMLNode& XML, DYNAMO::SimData* tmp):
  CICapture(tmp,NULL) //A temporary value!
{
  operator<<(XML);
}

void 
CISoftCore::operator<<(const XMLNode& XML)
{
  if (strcmp(XML.getAttribute("Type"),"SoftCore"))
    D_throw() << "Attempting to load SoftCore from non SoftCore entry";
  
  range.set_ptr(C2Range::loadClass(XML,Sim));
  
  try {
    diameter = Sim->Dynamics.units().unitLength() 
      * boost::lexical_cast<Iflt>(XML.getAttribute("Diameter"));
    
    wellDepth = boost::lexical_cast<Iflt>(XML.getAttribute("WellDepth"))
      * Sim->Dynamics.units().unitEnergy();
    
    d2 = diameter * diameter;
    
    intName = XML.getAttribute("Name");

    CICapture::loadCaptureMap(XML);   
  }
  catch (boost::bad_lexical_cast &)
    { D_throw() << "Failed a lexical cast in CISoftCore"; }
}

CInteraction* 
CISoftCore::Clone() const 
{ return new CISoftCore(*this); }

Iflt 
CISoftCore::hardCoreDiam() const 
{ return diameter; }

Iflt 
CISoftCore::maxIntDist() const 
{ return diameter; }

void 
CISoftCore::rescaleLengths(Iflt scale) 
{ 
  diameter += scale*diameter; 
  d2 = diameter*diameter;
}

void 
CISoftCore::initialise(size_t nID)
{
  ID = nID;
  CICapture::initCaptureMap();
}

bool 
CISoftCore::captureTest(const CParticle& p1, const CParticle& p2) const
{
  CVector<> rij = p1.getPosition() - p2.getPosition();
  Sim->Dynamics.BCs().setPBC(rij);
  
  return ((rij % rij) <= d2);
}

CIntEvent 
CISoftCore::getEvent(const CParticle &p1, 
		     const CParticle &p2) const 
{
#ifdef DYNAMO_DEBUG
  if (!Sim->Dynamics.Liouvillean().isUpToDate(p1))
    D_throw() << "Particle 1 is not up to date";
  
  if (!Sim->Dynamics.Liouvillean().isUpToDate(p2))
    D_throw() << "Particle 2 is not up to date";

  if (p1 == p2)
    D_throw() << "You shouldn't pass p1==p2 events to the interactions!";
#endif 

  CPDData colldat(*Sim, p1, p2);
    
  if (isCaptured(p1, p2)) 
    {
      if (Sim->Dynamics.Liouvillean().SphereSphereOutRoot(colldat, d2))
	return CIntEvent(p1, p2, colldat.dt, WELL_OUT, *this);
    }
  else if (Sim->Dynamics.Liouvillean().SphereSphereInRoot(colldat, d2)) 
    {
#ifdef DYNAMO_OverlapTesting
      if (Sim->Dynamics.Liouvillean().sphereOverlap(colldat,d2))
	D_throw() << "Overlapping cores (but not registered as captured) particles found in soft core" 
		  << "\nparticle1 " << p1.getID() << ", particle2 " 
		  << p2.getID() << "\nOverlap = " 
		  << (sqrt(colldat.r2) - sqrt(d2)) / Sim->Dynamics.units().unitLength();
#endif

      return CIntEvent(p1, p2, colldat.dt, WELL_IN, *this);
    }

  return CIntEvent(p1, p2, HUGE_VAL, NONE, *this);
}

void
CISoftCore::runEvent(const CParticle& p1, const CParticle& p2) const
{
  Sim->Dynamics.Liouvillean().updateParticlePair(p1, p2);
  CIntEvent iEvent = getEvent(p1, p2);
  
  if (iEvent.getType() == NONE)
    {
      I_cerr() << "A glancing or tenuous collision may have become invalid due"
	"\nto free streaming inaccuracies"
	"\nOtherwise the simulation has run out of events!"
	"\nThis occured when confirming the event with the scheduler"
	"\nIgnoring this NONE event below\n"
	       << iEvent.stringData(Sim);
      
      //Now we're past the event, update the scheduler and plugins
      Sim->ptrScheduler->fullUpdate(p1, p2);
      return;
    }

  ++Sim->lNColl;
  
#ifdef DYNAMO_DEBUG 
  if (isnan(iEvent.getdt()))
    D_throw() << "A NAN Interaction collision time has been found"
	      << iEvent.stringData(Sim);
  
  if (iEvent.getdt() == HUGE_VAL)
    D_throw() << "An infinite Interaction (not marked as NONE) collision time has been found\n"
	      << iEvent.stringData(Sim);
#endif
  
  //Debug section
#ifdef DYNAMO_CollDebug
  if (p2.getID() < p2.getID())
    std::cerr << "\nsysdt " << iEvent.getdt() + dSysTime
	      << "  ID1 " << p1.getID() 
	      << "  ID2 " << p2.getID()
	      << "  dt " << iEvent.getdt()
	      << "  Type " << CIntEvent::getCollEnumName(iEvent.getType());
  else
    std::cerr << "\nsysdt " << iEvent.getdt() + dSysTime
	      << "  ID1 " << p2().getID() 
	      << "  ID2 " << p1().getID()
	      << "  dt " << iEvent.getdt()
	      << "  Type " << CIntEvent::getCollEnumName(iEvent.getType());
#endif
  
  Sim->dSysTime += iEvent.getdt();
    
  Sim->ptrScheduler->stream(iEvent.getdt());
  
  //dynamics must be updated first
  Sim->Dynamics.stream(iEvent.getdt());

  switch (iEvent.getType())
    {
    case WELL_IN:
      {
	C2ParticleData retVal(Sim->Dynamics.Liouvillean()
			      .SphereWellEvent(iEvent, wellDepth, d2));
	
	if (retVal.getType() != BOUNCE)
	  addToCaptureMap(p1, p2);      
	
	//Now we're past the event, update the scheduler and plugins
	Sim->signalParticleUpdate(retVal);
	Sim->ptrScheduler->fullUpdate(p1, p2);
	
	BOOST_FOREACH(smrtPlugPtr<COutputPlugin> & Ptr, Sim->outputPlugins)
	  Ptr->eventUpdate(iEvent, retVal);


	break;
      }
    case WELL_OUT:
      {
	C2ParticleData retVal(Sim->Dynamics.Liouvillean()
			      .SphereWellEvent(iEvent, -wellDepth, d2));
	
	if (retVal.getType() != BOUNCE)
	  removeFromCaptureMap(p1, p2);      
	
	Sim->signalParticleUpdate(retVal);

	//Now we're past the event, update the scheduler and plugins
	Sim->ptrScheduler->fullUpdate(p1, p2);
	
	BOOST_FOREACH(smrtPlugPtr<COutputPlugin>& Ptr, 
		      Sim->outputPlugins)
	  Ptr->eventUpdate(iEvent, retVal);

	break;
      }
    default:
      D_throw() << "Unknown collision type";
    }
}

void
CISoftCore::checkOverlaps(const CParticle& part1, const CParticle& part2) const
{
  CVector<> rij = part1.getPosition() - part2.getPosition();
  Sim->Dynamics.BCs().setPBC(rij);
  Iflt r2 = rij.square();

  if (isCaptured(part1, part2))
    {
      if (r2 > d2)
	I_cerr() << std::setprecision(std::numeric_limits<float>::digits10)
		 << "Possible escaped captured pair in diagnostics\n ID1=" << part1.getID() 
		 << ", ID2=" << part2.getID() << "\nR_ij^2=" 
		 << r2 / pow(Sim->Dynamics.units().unitLength(),2)
		 << "\nd^2=" 
		 << d2 / pow(Sim->Dynamics.units().unitLength(),2);
    }
  else 
    if (r2 < d2)
      I_cerr() << std::setprecision(std::numeric_limits<float>::digits10)
	       << "Possible missed captured pair in diagnostics\n ID1=" << part1.getID() 
	       << ", ID2=" << part2.getID() << "\nR_ij^2=" 
	       << r2 / pow(Sim->Dynamics.units().unitLength(),2)
	       << "\nd^2=" 
	       << d2 / pow(Sim->Dynamics.units().unitLength(),2);
}
  
void 
CISoftCore::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "SoftCore"
      << xmlw::attr("Diameter") 
      << diameter / Sim->Dynamics.units().unitLength() 
      << xmlw::attr("WellDepth") 
      << wellDepth / Sim->Dynamics.units().unitEnergy()
      << xmlw::attr("Name") << intName
      << range;
  
  CICapture::outputCaptureMap(XML);  
}
