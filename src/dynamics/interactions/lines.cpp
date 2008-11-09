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

#include "lines.hpp"
#include <boost/lexical_cast.hpp>
#include <cmath>
#include <iomanip>
#include "../../base/is_exception.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../../dynamics/interactions/intEvent.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../liouvillean/OrientationL.hpp"
#include "../units/units.hpp"
#include "../../base/is_simdata.hpp"
#include "../2particleEventData.hpp"
#include "../BC/BC.hpp"

CILines::CILines(const DYNAMO::SimData* tmp, Iflt nd, 
		 Iflt ne, C2Range* nR):
  CICapture(tmp, nR),
  length(nd), l2(nd*nd), e(ne) {}

CILines::CILines(const XMLNode& XML, const DYNAMO::SimData* tmp):
  CICapture(tmp, NULL)
{
  operator<<(XML);
}

void 
CILines::initialise(size_t nID)
{
  if (dynamic_cast<const CLNOrientation*>(&(Sim->Dynamics.Liouvillean()))
      == NULL)
    D_throw() << "Interaction requires an orientation capable Liouvillean.";

  ID = nID; 

  CICapture::initCaptureMap();
}

void 
CILines::operator<<(const XMLNode& XML)
{ 
  if (strcmp(XML.getAttribute("Type"),"Lines"))
    D_throw() << "Attempting to load Lines from non Lines entry";
  
  range.set_ptr(C2Range::loadClass(XML,Sim));
  
  try 
    {
      length = Sim->Dynamics.units().unitLength() * 
	boost::lexical_cast<Iflt>(XML.getAttribute("Length"));
      
      l2 = length * length;

      e = boost::lexical_cast<Iflt>(XML.getAttribute("Elasticity"));
      
      intName = XML.getAttribute("Name");

      CICapture::loadCaptureMap(XML);   
    }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in CILines";
    }
}

Iflt 
CILines::maxIntDist() const 
{ return length; }

Iflt 
CILines::hardCoreDiam() const 
{ return 0.0; }

void 
CILines::rescaleLengths(Iflt scale) 
{ 
  length += scale * length;
  
  l2 = length * length;
}

CInteraction* 
CILines::Clone() const 
{ return new CILines(*this); }
  
CIntEvent 
CILines::getCollision(const CParticle &p1, const CParticle &p2) const 
{ 
  //We use wells to mark when to do the line test
  Sim->Dynamics.Liouvillean().updateParticlePair(p1, p2);
  CPDData colldat(*Sim, p1, p2);
  
  if (isCaptured(p1, p2)) 
    {
      //Run this to determine when the spheres no longer intersect
      if (!Sim->Dynamics.Liouvillean().SphereSphereOutRoot(colldat, l2))
	//The spheres never stop intersecting!
	D_throw() << "Error, spheres always remain overlapping? "
	  "\ncompressed lines? Not implemented yet";

      //colldat.dt has the upper limit of the line collision time
      //Lower limit is right now
      //Test for a line collision
      if (Sim->Dynamics.Liouvillean().getLineLineCollision
	  (colldat, length, p1, p2, colldat.dt))
	return CIntEvent(p1, p2, colldat.dt, CORE, *this);

      return CIntEvent(p1, p2, colldat.dt, WELL_OUT, *this);
    }
  else if (Sim->Dynamics.Liouvillean().SphereSphereInRoot(colldat, l2)) 
    return CIntEvent(p1, p2, colldat.dt, WELL_IN, *this);
  
  return CIntEvent(p1, p2, HUGE_VAL, NONE, *this);
}

C2ParticleData 
CILines::runCollision(const CIntEvent &event) const
{ 
  switch (event.getType())
    {
    case CORE:
      {
	//We have a line interaction! Run it
	return Sim->Dynamics.Liouvillean().runLineLineCollision(event);
      }
    case WELL_IN:
      {	
	addToCaptureMap(event.getParticle1(), event.getParticle2());      
	
	C2ParticleData retval(event.getParticle1(), event.getParticle2(), 
			      Sim->Dynamics.getSpecies(event.getParticle1()),
			      Sim->Dynamics.getSpecies(event.getParticle2()),
			      VIRTUAL);
	
	retval.dP = CVector<>(0.0);
	retval.deltake = CVector<>(0.0);	

	return retval;
      }
    case WELL_OUT:
      {
	removeFromCaptureMap(event.getParticle1(), event.getParticle2());      
	
	C2ParticleData retval(event.getParticle1(), event.getParticle2(), 
			      Sim->Dynamics.getSpecies(event.getParticle1()),
			      Sim->Dynamics.getSpecies(event.getParticle2()),
			      VIRTUAL);
	
	retval.dP = CVector<>(0.0);
	retval.deltake = CVector<>(0.0);	

	return retval;
      }
    default:
      D_throw() << "Unknown collision type";
    }
}
   
void 
CILines::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "Lines"
      << xmlw::attr("Length") << length / Sim->Dynamics.units().unitLength()
      << xmlw::attr("Elasticity") << e
      << xmlw::attr("Name") << intName
      << range;

  CICapture::outputCaptureMap(XML);
}

bool 
CILines::captureTest(const CParticle& p1, const CParticle& p2) const
{
  CVector<> rij = p1.getPosition() - p2.getPosition();
  Sim->Dynamics.BCs().setPBC(rij);
  
  return rij % rij <= l2;
}

void
CILines::checkOverlaps(const CParticle& part1, const CParticle& part2) const
{}
