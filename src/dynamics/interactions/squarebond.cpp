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

#include "squarebond.hpp"
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
#include "../../schedulers/scheduler.hpp"
#include "../NparticleEventData.hpp"
#include <iomanip>

CISquareBond::CISquareBond(DYNAMO::SimData* tmp, Iflt nd, Iflt nl, C2Range* nR):
  CInteraction(tmp,nR),
  diameter(nd),d2(nd*nd),lambda(nl),ld2(nd*nd*nl*nl) 
{}
  
CISquareBond::CISquareBond(const XMLNode& XML, DYNAMO::SimData* tmp):
  CInteraction(tmp,NULL) //A temporary value!
{ operator<<(XML); }
	    
void 
CISquareBond::operator<<(const XMLNode& XML)
{
if (strcmp(XML.getAttribute("Type"),"SquareBond"))
  D_throw() << "Attempting to load SquareBond from non SquareBond entry";
      
    range.set_ptr(C2Range::loadClass(XML,Sim));
      
      try {
    diameter = boost::lexical_cast<Iflt>(XML.getAttribute("Diameter"))
      * Sim->Dynamics.units().unitLength();

    lambda = boost::lexical_cast<Iflt>(XML.getAttribute("Lambda"));

    intName = XML.getAttribute("Name");
    
    d2 = diameter * diameter;

    ld2 = d2 * lambda * lambda;
  }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in CISquareWell";
    }
}

CInteraction* 
CISquareBond::Clone() const 
{ return new CISquareBond(*this); }

Iflt 
CISquareBond::getCaptureEnergy() const 
{ return 0.0; }

Iflt 
CISquareBond::hardCoreDiam() const 
{ return diameter; }

Iflt 
CISquareBond::maxIntDist() const 
{ return diameter*lambda; }

void 
CISquareBond::rescaleLengths(Iflt scale) 
{ 
  diameter += scale*diameter;
  d2 = diameter*diameter;
  ld2 = diameter*diameter*lambda*lambda;
}

void 
CISquareBond::initialise(size_t nID)
{
  ID = nID;
}

bool 
CISquareBond::captureTest(const CParticle& p1, const CParticle& p2) const
{
  Vector  rij = p1.getPosition() - p2.getPosition();
  Sim->Dynamics.BCs().applyBC(rij);
  
  if (((rij | rij) <= ld2) && ((rij | rij) >= d2))
    return true;
  
  return false;
}

void
CISquareBond::checkOverlaps(const CParticle& part1, const CParticle& part2) const
{
  Vector  rij = part1.getPosition() - part2.getPosition();
  Sim->Dynamics.BCs().applyBC(rij);
  Iflt r2 = rij.nrm2();

  if (r2 < d2)
    I_cerr() << std::setprecision(std::numeric_limits<float>::digits10)
	     << "Possible bonded overlap occured in diagnostics\n ID1=" << part1.getID() 
	     << ", ID2=" << part2.getID() << "\nR_ij^2=" 
	     << r2 / pow(Sim->Dynamics.units().unitLength(),2)
	     << "\nd^2=" 
	     << d2 / pow(Sim->Dynamics.units().unitLength(),2);
  
  if (r2 > ld2)
    I_cerr() << std::setprecision(std::numeric_limits<float>::digits10)
	     << "Possible escaped bonded pair in diagnostics\n ID1=" << part1.getID() 
	     << ", ID2=" << part2.getID() << "\nR_ij^2=" 
	     << r2 / pow(Sim->Dynamics.units().unitLength(),2)
	     << "\n(lambda * d)^2=" 
	     << ld2 / pow(Sim->Dynamics.units().unitLength(),2);
}

CIntEvent 
CISquareBond::getEvent(const CParticle &p1, 
		       const CParticle &p2) const 
{    
#ifdef DYNAMO_DEBUG
  if (!Sim->Dynamics.getLiouvillean().isUpToDate(p1))
    D_throw() << "Particle 1 is not up to date";
  
  if (!Sim->Dynamics.getLiouvillean().isUpToDate(p2))
    D_throw() << "Particle 2 is not up to date";

  if (p1 == p2)
    D_throw() << "You shouldn't pass p1==p2 events to the interactions!";
#endif 

  CPDData colldat(*Sim, p1, p2);

  if (Sim->Dynamics.getLiouvillean().SphereSphereInRoot(colldat, d2))
    {
#ifdef DYNAMO_OverlapTesting
      if (Sim->Dynamics.getLiouvillean().sphereOverlap(colldat,d2))
	D_throw() << "Overlapping particles found" 
		  << ", particle1 " << p1.getID() 
		  << ", particle2 " 
		  << p2.getID() << "\nOverlap = " << (sqrt(colldat.r2) - sqrt(d2))/Sim->Dynamics.units().unitLength();
#endif      
      return CIntEvent(p1, p2, colldat.dt, CORE, *this);
    }
  else
    if (Sim->Dynamics.getLiouvillean().SphereSphereOutRoot(colldat, ld2))
      {
	return CIntEvent(p1, p2, colldat.dt, BOUNCE, *this); 
      }
  
  return CIntEvent(p1, p2, HUGE_VAL, NONE, *this);
}

void
CISquareBond::runEvent(const CParticle& p1, const CParticle& p2, 
		       const CIntEvent& iEvent) const
{
  ++Sim->lNColl;

#ifdef DYNAMO_DEBUG
  if ((iEvent.getType() != BOUNCE) && (iEvent.getType() != CORE))
    D_throw() << "Unknown type found";
#endif

  C2ParticleData EDat(Sim->Dynamics.getLiouvillean().SmoothSpheresColl
		      (iEvent, 1.0, d2, iEvent.getType()));

  Sim->signalParticleUpdate(EDat);
    
  //Now we're past the event, update the scheduler and plugins
  Sim->ptrScheduler->fullUpdate(p1, p2);
  
  BOOST_FOREACH(smrtPlugPtr<COutputPlugin> & Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(iEvent,EDat);

}
    
void 
CISquareBond::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "SquareBond"
      << xmlw::attr("Diameter") 
      << diameter / Sim->Dynamics.units().unitLength()
      << xmlw::attr("Lambda") << lambda
      << xmlw::attr("Name") << intName
      << range;
}

void 
CISquareBond::write_povray_info(std::ostream& os) const
{  
  BOOST_FOREACH(const CParticle& p1, Sim->vParticleList)
    BOOST_FOREACH(const CParticle& p2, Sim->vParticleList)
    if (range->isInRange(p1,p2) && (p1 != p2))
      {
	Vector  pos1(p1.getPosition()), pos2(p2.getPosition());
	Sim->Dynamics.BCs().applyBC(pos1);
	Sim->Dynamics.BCs().applyBC(pos2);
	
	if ((pos1-pos2).nrm() > 0.5) continue;

	Sim->Dynamics.BCs().applyBC(pos1);
	Sim->Dynamics.BCs().applyBC(pos2);
	
	os << "cylinder {\n <"
	   << pos1[0];
	
	for (size_t iDim(1); iDim < NDIM; ++iDim)
	  os << "," << pos1[iDim];
	
	os << ">, <" << pos2[0];
	
	for (size_t iDim(1); iDim < NDIM; ++iDim)
	  os << "," << pos2[iDim];
	
	os << ">, " << 0.1 * diameter << " pigment{color Green}}\n";
      }

}
