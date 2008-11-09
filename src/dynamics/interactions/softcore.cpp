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

CISoftCore::CISoftCore(const DYNAMO::SimData* tmp, Iflt nd, Iflt nWD, 
		       C2Range* nR):
  CICapture(tmp,nR),
  diameter(nd),d2(nd*nd),wellDepth(nWD) 
{}

CISoftCore::CISoftCore(const XMLNode& XML, const DYNAMO::SimData* tmp):
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
CISoftCore::getCollision(const CParticle &p1, 
			 const CParticle &p2) const 
{
#ifdef DYNAMO_DEBUG
  if (p1 == p2)
    D_throw() << "You shouldn't pass p1==p2 events to the interactions!";
#endif 

  Sim->Dynamics.Liouvillean().updateParticlePair(p1, p2);
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

C2ParticleData 
CISoftCore::runCollision(const CIntEvent &event) const
{  
  switch (event.getType())
    {
    case WELL_IN:
      {
	C2ParticleData retVal(Sim->Dynamics.Liouvillean()
			      .SphereWellEvent(event, wellDepth, d2));
	
	if (retVal.getType() != BOUNCE)
	  addToCaptureMap(event.getParticle1(), event.getParticle2());      
	
	return retVal;
      }
    case WELL_OUT:
      {
	C2ParticleData retVal(Sim->Dynamics.Liouvillean()
			      .SphereWellEvent(event, -wellDepth, d2));
	
	if (retVal.getType() != BOUNCE)
	  removeFromCaptureMap(event.getParticle1(), event.getParticle2());      
	
	return retVal;
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
  double r2 = rij.square();

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
