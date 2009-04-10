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

#include "squarewell.hpp"
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

CISquareWell::CISquareWell(DYNAMO::SimData* tmp, Iflt nd, Iflt nl, 
			   Iflt nWD, 
			   Iflt ne, C2Range* nR):
  CICapture(tmp,nR),
  diameter(nd), d2(nd*nd), lambda(nl), 
  ld2(nd*nd*nl*nl), wellDepth(nWD),
  e(ne) {}

CISquareWell::CISquareWell(const XMLNode& XML, DYNAMO::SimData* tmp):
  CICapture(tmp, NULL) //A temporary value!
{
  operator<<(XML);
}

void 
CISquareWell::operator<<(const XMLNode& XML)
{
  if (strcmp(XML.getAttribute("Type"),"SquareWell"))
    D_throw() << "Attempting to load SquareWell from non SquareWell entry";
  
  range.set_ptr(C2Range::loadClass(XML,Sim));
  
  try {
    diameter = Sim->Dynamics.units().unitLength() 
      * boost::lexical_cast<Iflt>(XML.getAttribute("Diameter"));
    
    e = boost::lexical_cast<Iflt>(XML.getAttribute("Elasticity"));
    
    wellDepth = boost::lexical_cast<Iflt>(XML.getAttribute("WellDepth"))
      * Sim->Dynamics.units().unitEnergy();
    
    lambda = boost::lexical_cast<Iflt>(XML.getAttribute("Lambda"));
    
    d2 = diameter * diameter;
    
    ld2 = d2 * lambda * lambda;
    
    intName = XML.getAttribute("Name");

    CICapture::loadCaptureMap(XML);   
  }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in CISquareWell";
    }
}

CInteraction* 
CISquareWell::Clone() const 
{ return new CISquareWell(*this); }

Iflt 
CISquareWell::hardCoreDiam() const 
{ return diameter; }

Iflt 
CISquareWell::maxIntDist() const 
{ return diameter*lambda; }

void 
CISquareWell::rescaleLengths(Iflt scale) 
{ 
  diameter += scale*diameter; 

  d2 = diameter*diameter;

  ld2 = diameter*diameter*lambda*lambda;
}

void 
CISquareWell::initialise(size_t nID)
{
  ID = nID;
  CICapture::initCaptureMap();
}

bool 
CISquareWell::captureTest(const CParticle& p1, const CParticle& p2) const
{
  CVector<> rij = p1.getPosition() - p2.getPosition();
  Sim->Dynamics.BCs().setPBC(rij);
  
  return (((rij | rij) <= ld2) && ((rij | rij) >= d2));
}

CIntEvent
CISquareWell::getEvent(const CParticle &p1, 
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
      if (Sim->Dynamics.Liouvillean().SphereSphereInRoot(colldat, d2)) 
	{
#ifdef DYNAMO_OverlapTesting
	  //Check that there is no overlap 
	  if (Sim->Dynamics.Liouvillean().sphereOverlap(colldat, d2))
	    D_throw() << "Overlapping particles found" 
		      << ", particle1 " << p1.getID() 
		      << ", particle2 " 
		      << p2.getID() << "\nOverlap = " << (sqrt(colldat.r2) - sqrt(d2))/Sim->Dynamics.units().unitLength();
#endif
	  
	  return CIntEvent(p1, p2, colldat.dt, CORE, *this);
	}
      else
	if (Sim->Dynamics.Liouvillean().SphereSphereOutRoot(colldat, ld2))
	  {  
	    return CIntEvent(p1, p2, colldat.dt, WELL_OUT, *this);
	  }
    }
  else if (Sim->Dynamics.Liouvillean().SphereSphereInRoot(colldat, ld2)) 
    {
#ifdef DYNAMO_OverlapTesting
      if (Sim->Dynamics.Liouvillean().sphereOverlap(colldat,ld2))
	{
	  if (Sim->Dynamics.Liouvillean().sphereOverlap(colldat,d2))
	    D_throw() << "Overlapping cores (but not registerd as captured) particles found in square well" 
		      << "\nparticle1 " << p1.getID() << ", particle2 " 
		      << p2.getID() << "\nOverlap = " 
		      << (sqrt(colldat.r2) - sqrt(d2)) 
	      / Sim->Dynamics.units().unitLength();
	  else
	    D_throw() << "Overlapping wells (but not registerd as captured) particles found" 
		      << "\nparticle1 " << p1.getID() << ", particle2 " 
		      << p2.getID() << "\nOverlap = " 
		      << (sqrt(colldat.r2) - sqrt(ld2)) 
	      / Sim->Dynamics.units().unitLength();
	  
	}
#endif

      return CIntEvent(p1, p2, colldat.dt, WELL_IN, *this);
    }

  return CIntEvent(p1, p2, HUGE_VAL, NONE, *this);
}

void
CISquareWell::runEvent(const CParticle& p1, 
		       const CParticle& p2,
		       const CIntEvent& iEvent) const
{
  ++Sim->lNColl;

  switch (iEvent.getType())
    {
    case CORE:
      {
	C2ParticleData retVal(Sim->Dynamics.Liouvillean().SmoothSpheresColl(iEvent, e, d2, CORE));
	Sim->signalParticleUpdate(retVal);
	
	Sim->ptrScheduler->fullUpdate(p1, p2);
	
	BOOST_FOREACH(smrtPlugPtr<COutputPlugin> & Ptr, Sim->outputPlugins)
	  Ptr->eventUpdate(iEvent, retVal);

	break;
      }
    case WELL_IN:
      {
	C2ParticleData retVal(Sim->Dynamics.Liouvillean()
			      .SphereWellEvent(iEvent, wellDepth, ld2));
	
	if (retVal.getType() != BOUNCE)
	  addToCaptureMap(p1, p2);      
	
	Sim->ptrScheduler->fullUpdate(p1, p2);
	Sim->signalParticleUpdate(retVal);
	
	BOOST_FOREACH(smrtPlugPtr<COutputPlugin> & Ptr, Sim->outputPlugins)
	  Ptr->eventUpdate(iEvent, retVal);


	break;
      }
    case WELL_OUT:
      {
	C2ParticleData retVal(Sim->Dynamics.Liouvillean()
			      .SphereWellEvent(iEvent, -wellDepth, ld2));
	
	if (retVal.getType() != BOUNCE)
	  removeFromCaptureMap(p1, p2);      

	Sim->signalParticleUpdate(retVal);

	Sim->ptrScheduler->fullUpdate(p1, p2);
	
	BOOST_FOREACH(smrtPlugPtr<COutputPlugin> & Ptr, Sim->outputPlugins)
	  Ptr->eventUpdate(iEvent, retVal);
	break;
      }
    default:
      D_throw() << "Unknown collision type";
    } 
}

void
CISquareWell::checkOverlaps(const CParticle& part1, const CParticle& part2) const
{
  CVector<> rij = part1.getPosition() - part2.getPosition();
  Sim->Dynamics.BCs().setPBC(rij);
  Iflt r2 = rij.square();

  if (isCaptured(part1, part2))
    {
      if (r2 < d2)
	I_cerr() << std::setprecision(std::numeric_limits<float>::digits10)
		 << "Possible captured overlap occured in diagnostics\n ID1=" << part1.getID() 
		 << ", ID2=" << part2.getID() << "\nR_ij^2=" 
		 << r2 / pow(Sim->Dynamics.units().unitLength(),2)
		 << "\nd^2=" 
		 << d2 / pow(Sim->Dynamics.units().unitLength(),2);

      if (r2 > lambda * lambda * d2)
	I_cerr() << std::setprecision(std::numeric_limits<float>::digits10)
		 << "Possible escaped captured pair in diagnostics\n ID1=" << part1.getID() 
		 << ", ID2=" << part2.getID() << "\nR_ij^2=" 
		 << r2 / pow(Sim->Dynamics.units().unitLength(),2)
		 << "\n(lambda * d)^2=" 
		 << ld2 / pow(Sim->Dynamics.units().unitLength(),2);
    }
  else 
    if (r2 < ld2)
      I_cerr() << std::setprecision(std::numeric_limits<float>::digits10)
	       << "Possible missed captured pair in diagnostics\n ID1=" << part1.getID() 
	       << ", ID2=" << part2.getID() << "\nR_ij^2=" 
	       << r2 / pow(Sim->Dynamics.units().unitLength(),2)
	       << "\n(lambda * d)^2=" 
	       << ld2 / pow(Sim->Dynamics.units().unitLength(),2);
}
  
void 
CISquareWell::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "SquareWell"
      << xmlw::attr("Diameter") 
      << diameter / Sim->Dynamics.units().unitLength() 
      << xmlw::attr("Elasticity") << e
      << xmlw::attr("Lambda") << lambda
      << xmlw::attr("WellDepth") 
      << wellDepth / Sim->Dynamics.units().unitEnergy()
      << xmlw::attr("Name") << intName
      << range;
  
  CICapture::outputCaptureMap(XML);  
}

void 
CISquareWell::write_povray_desc(const DYNAMO::RGB& rgb, 
				const size_t& specID, 
				std::ostream& os) const
{
  os << "#declare intrep" << ID << "center = " 
     << "sphere {\n <0,0,0> " 
     << diameter / 2.0
     << "\n texture { pigment { color rgb<" << rgb.R << "," << rgb.G 
     << "," << rgb.B << "> }}\nfinish { phong 0.9 phong_size 60 }\n}\n"
     << "#declare intrep" << ID << "well = sphere {\n <0,0,0> " << diameter * lambda * 0.5
     << "\n texture { pigment { color rgbt <1,1,1,0.9> }}\n}\n";

  BOOST_FOREACH(const size_t& part, *(Sim->Dynamics.getSpecies()[specID].getRange()))
    {
      CVector<> pos(Sim->vParticleList[part].getPosition());
      Sim->Dynamics.BCs().setPBC(pos);
      
      os << "object {\n intrep" << ID << "center\n translate < "
	 << pos[0] << ", " << pos[1] << ", " << pos[2] << ">\n}\n";
    }
  
  os << "merge {\n";
  BOOST_FOREACH(const size_t& part, *(Sim->Dynamics.getSpecies()[specID].getRange()))
    {
      CVector<> pos(Sim->vParticleList[part].getPosition());
      Sim->Dynamics.BCs().setPBC(pos);

      os << "object {\n intrep" << ID << "well\n translate < "
	 << pos[0] << ", " << pos[1] << ", " << pos[2] << ">\n}\n";
    }
  
  os << "}\n";

}
