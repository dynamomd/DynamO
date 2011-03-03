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

#include "squarewell.hpp"
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
#include "../../schedulers/scheduler.hpp"
#include "../NparticleEventData.hpp"
#include <iomanip>

ISquareWell::ISquareWell(DYNAMO::SimData* tmp, double nd, double nl, 
			   double nWD, 
			   double ne, C2Range* nR):
  ISingleCapture(tmp,nR),
  diameter(nd), d2(nd*nd), lambda(nl), 
  ld2(nd*nd*nl*nl), wellDepth(nWD),
  e(ne) {}

ISquareWell::ISquareWell(const XMLNode& XML, DYNAMO::SimData* tmp):
  ISingleCapture(tmp, NULL) //A temporary value!
{
  operator<<(XML);
}

void 
ISquareWell::operator<<(const XMLNode& XML)
{
  if (strcmp(XML.getAttribute("Type"),"SquareWell"))
    M_throw() << "Attempting to load SquareWell from non SquareWell entry";
  
  range.set_ptr(C2Range::loadClass(XML,Sim));
  
  try {
    diameter = Sim->dynamics.units().unitLength() 
      * boost::lexical_cast<double>(XML.getAttribute("Diameter"));
    
    e = boost::lexical_cast<double>(XML.getAttribute("Elasticity"));
    
    wellDepth = boost::lexical_cast<double>(XML.getAttribute("WellDepth"))
      * Sim->dynamics.units().unitEnergy();
    
    lambda = boost::lexical_cast<double>(XML.getAttribute("Lambda"));
    
    d2 = diameter * diameter;
    
    ld2 = d2 * lambda * lambda;
    
    intName = XML.getAttribute("Name");

    ISingleCapture::loadCaptureMap(XML);   
  }
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in CISquareWell";
    }
}

Interaction* 
ISquareWell::Clone() const 
{ return new ISquareWell(*this); }

double 
ISquareWell::hardCoreDiam() const 
{ return diameter; }

double 
ISquareWell::maxIntDist() const 
{ return diameter * lambda; }

void 
ISquareWell::rescaleLengths(double scale) 
{ 
  diameter += scale*diameter; 

  d2 = diameter*diameter;

  ld2 = diameter*diameter*lambda*lambda;
}

void 
ISquareWell::initialise(size_t nID)
{
  ID = nID;
  ISingleCapture::initCaptureMap();
}

bool 
ISquareWell::captureTest(const Particle& p1, const Particle& p2) const
{
  Vector  rij = p1.getPosition() - p2.getPosition();
  Sim->dynamics.BCs().applyBC(rij);

#ifdef DYNAMO_DEBUG
  if ((rij | rij) < d2)
    I_cerr() << "Warning! Two particles might be overlapping"
	     << "\nrij^2 = " << (rij | rij)
	     << "\nd^2 = " << d2;
#endif

  return (rij | rij) <= ld2;
}

IntEvent
ISquareWell::getEvent(const Particle &p1, 
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
      if (Sim->dynamics.getLiouvillean()
	  .SphereSphereInRoot(colldat, d2,
			      p1.testState(Particle::DYNAMIC), p2.testState(Particle::DYNAMIC))) 
	{
#ifdef DYNAMO_OverlapTesting
	  //Check that there is no overlap 
	  if (Sim->dynamics.getLiouvillean().sphereOverlap(colldat, d2))
	    M_throw() << "Overlapping particles found" 
		      << ", particle1 " << p1.getID() 
		      << ", particle2 " 
		      << p2.getID() << "\nOverlap = " << (sqrt(colldat.r2) - sqrt(d2))/Sim->dynamics.units().unitLength();
#endif
	  
	  return IntEvent(p1, p2, colldat.dt, CORE, *this);
	}
      else
	if (Sim->dynamics.getLiouvillean()
	    .SphereSphereOutRoot(colldat, ld2,
				 p1.testState(Particle::DYNAMIC), p2.testState(Particle::DYNAMIC)))
	  {  
	    return IntEvent(p1, p2, colldat.dt, WELL_OUT, *this);
	  }
    }
  else if (Sim->dynamics.getLiouvillean()
	   .SphereSphereInRoot(colldat, ld2, 
			       p1.testState(Particle::DYNAMIC), p2.testState(Particle::DYNAMIC))) 
    {
#ifdef DYNAMO_OverlapTesting
      if (Sim->dynamics.getLiouvillean().sphereOverlap(colldat,ld2))
	{
	  if (Sim->dynamics.getLiouvillean().sphereOverlap(colldat,d2))
	    M_throw() << "Overlapping cores (but not registerd as captured) particles found in square well" 
		      << "\nparticle1 " << p1.getID() << ", particle2 " 
		      << p2.getID() << "\nOverlap = " 
		      << (sqrt(colldat.r2) - sqrt(d2)) 
	      / Sim->dynamics.units().unitLength();
	  else
	    M_throw() << "Overlapping wells (but not registerd as captured) particles found" 
		      << "\nparticle1 " << p1.getID() << ", particle2 " 
		      << p2.getID() << "\nOverlap = " 
		      << (sqrt(colldat.r2) - sqrt(ld2)) 
	      / Sim->dynamics.units().unitLength();
	  
	}
#endif

      return IntEvent(p1, p2, colldat.dt, WELL_IN, *this);
    }

  return IntEvent(p1, p2, HUGE_VAL, NONE, *this);
}

void
ISquareWell::runEvent(const Particle& p1, 
		       const Particle& p2,
		       const IntEvent& iEvent) const
{
  ++Sim->eventCount;

  switch (iEvent.getType())
    {
    case CORE:
      {
	PairEventData retVal(Sim->dynamics.getLiouvillean().SmoothSpheresColl(iEvent, e, d2, CORE));
	Sim->signalParticleUpdate(retVal);
	
	Sim->ptrScheduler->fullUpdate(p1, p2);
	
	BOOST_FOREACH(magnet::ClonePtr<OutputPlugin> & Ptr, Sim->outputPlugins)
	  Ptr->eventUpdate(iEvent, retVal);

	break;
      }
    case WELL_IN:
      {
	PairEventData retVal(Sim->dynamics.getLiouvillean()
			      .SphereWellEvent(iEvent, wellDepth, ld2));
	
	if (retVal.getType() != BOUNCE)
	  addToCaptureMap(p1, p2);      
	
	Sim->ptrScheduler->fullUpdate(p1, p2);
	Sim->signalParticleUpdate(retVal);
	
	BOOST_FOREACH(magnet::ClonePtr<OutputPlugin> & Ptr, Sim->outputPlugins)
	  Ptr->eventUpdate(iEvent, retVal);


	break;
      }
    case WELL_OUT:
      {
	PairEventData retVal(Sim->dynamics.getLiouvillean()
			      .SphereWellEvent(iEvent, -wellDepth, ld2));
	
	if (retVal.getType() != BOUNCE)
	  removeFromCaptureMap(p1, p2);      

	Sim->signalParticleUpdate(retVal);

	Sim->ptrScheduler->fullUpdate(p1, p2);
	
	BOOST_FOREACH(magnet::ClonePtr<OutputPlugin> & Ptr, Sim->outputPlugins)
	  Ptr->eventUpdate(iEvent, retVal);
	break;
      }
    default:
      M_throw() << "Unknown collision type";
    } 
}

void
ISquareWell::checkOverlaps(const Particle& part1, const Particle& part2) const
{
  Vector  rij = part1.getPosition() - part2.getPosition();
  Sim->dynamics.BCs().applyBC(rij);
  double r2 = rij.nrm2();

  if (isCaptured(part1, part2))
    {
      if (r2 < d2)
	I_cerr() << std::setprecision(std::numeric_limits<float>::digits10)
		 << "Possible captured overlap occured in diagnostics\n ID1=" << part1.getID() 
		 << ", ID2=" << part2.getID() << "\nR_ij^2=" 
		 << r2 / pow(Sim->dynamics.units().unitLength(),2)
		 << "\nd^2=" 
		 << d2 / pow(Sim->dynamics.units().unitLength(),2);

      if (r2 > lambda * lambda * d2)
	I_cerr() << std::setprecision(std::numeric_limits<float>::digits10)
		 << "Possible escaped captured pair in diagnostics\n ID1=" << part1.getID() 
		 << ", ID2=" << part2.getID() << "\nR_ij^2=" 
		 << r2 / pow(Sim->dynamics.units().unitLength(),2)
		 << "\n(lambda * d)^2=" 
		 << ld2 / pow(Sim->dynamics.units().unitLength(),2);
    }
  else 
    if (r2 < ld2)
      {
	if (r2 < d2)
	  I_cerr() << std::setprecision(std::numeric_limits<float>::digits10)
		   << "Overlap error\n ID1=" << part1.getID() 
		   << ", ID2=" << part2.getID() << "\nR_ij^2=" 
		   << r2 / pow(Sim->dynamics.units().unitLength(),2)
		   << "\n(d)^2=" 
		   << d2 / pow(Sim->dynamics.units().unitLength(),2);
	else
	  I_cerr() << std::setprecision(std::numeric_limits<float>::digits10)
		   << "Possible missed captured pair in diagnostics\n ID1=" << part1.getID() 
		   << ", ID2=" << part2.getID() << "\nR_ij^2=" 
		   << r2 / pow(Sim->dynamics.units().unitLength(),2)
		   << "\n(lambda * d)^2=" 
		   << ld2 / pow(Sim->dynamics.units().unitLength(),2);
      }
}
  
void 
ISquareWell::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Type") << "SquareWell"
      << xml::attr("Diameter") 
      << diameter / Sim->dynamics.units().unitLength() 
      << xml::attr("Elasticity") << e
      << xml::attr("Lambda") << lambda
      << xml::attr("WellDepth") 
      << wellDepth / Sim->dynamics.units().unitEnergy()
      << xml::attr("Name") << intName
      << range;
  
  ISingleCapture::outputCaptureMap(XML);  
}

void 
ISquareWell::write_povray_desc(const DYNAMO::RGB& rgb, 
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

  BOOST_FOREACH(const size_t& part, *(Sim->dynamics.getSpecies()[specID]->getRange()))
    {
      Vector  pos(Sim->particleList[part].getPosition());
      Sim->dynamics.BCs().applyBC(pos);
      
      os << "object {\n intrep" << ID << "center\n translate < "
	 << pos[0] << ", " << pos[1] << ", " << pos[2] << ">\n}\n";
    }
  
  os << "merge {\n";
  BOOST_FOREACH(const size_t& part, *(Sim->dynamics.getSpecies()[specID]->getRange()))
    {
      Vector  pos(Sim->particleList[part].getPosition());
      Sim->dynamics.BCs().applyBC(pos);

      os << "object {\n intrep" << ID << "well\n translate < "
	 << pos[0] << ", " << pos[1] << ", " << pos[2] << ">\n}\n";
    }
  
  os << "}\n";

}
