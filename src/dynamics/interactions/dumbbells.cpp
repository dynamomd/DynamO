/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>
    Copyright (C) 2011  Sebastian Gonzalez <tsuresuregusa@gmail.com>

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

#include "dumbbells.hpp"
#include "../../dynamics/interactions/intEvent.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../liouvillean/OrientationL.hpp"
#include "../units/units.hpp"
#include "../../base/is_simdata.hpp"
#include "../2particleEventData.hpp"
#include "../BC/BC.hpp"
#include "../ranges/1range.hpp"
#include "../../schedulers/scheduler.hpp"
#include "../NparticleEventData.hpp"
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <cmath>
#include <iomanip>


IDumbbells::IDumbbells(DYNAMO::SimData* tmp, double nd, 
		       double ne, double ndiameter, C2Range* nR):
  ISingleCapture(tmp, nR),
  length(nd), e(ne), diameter(ndiameter) 
{}

IDumbbells::IDumbbells(const magnet::xml::Node& XML, DYNAMO::SimData* tmp):
  ISingleCapture(tmp, NULL)
{
  operator<<(XML);
}

void 
IDumbbells::initialise(size_t nID)
{
  if (dynamic_cast<const LNOrientation*>(&(Sim->dynamics.getLiouvillean()))
      == NULL)
    M_throw() << "Interaction requires an orientation capable Liouvillean.";
  
  ID = nID; 
  
  ISingleCapture::initCaptureMap();
}

void 
IDumbbells::operator<<(const magnet::xml::Node& XML)
{ 
  if (strcmp(XML.getAttribute("Type"),"Dumbbells"))
    M_throw() << "Attempting to load Dumbbells from non Dumbbells entry";
  
  range.set_ptr(C2Range::getClass(XML,Sim));
  
  try 
    {
      length = XML.getAttribute("Length").as<double>() * Sim->dynamics.units().unitLength();
      e = XML.getAttribute("Elasticity").as<double>();
      diameter = XML.getAttribute("Diameter").as<double>() * Sim->dynamics.units().unitLength();
      intName = XML.getAttribute("Name");
      ISingleCapture::loadCaptureMap(XML);   
    }
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in CIDumbbells";
    }
}

double 
IDumbbells::maxIntDist() const 
{ return length + diameter; }

double 
IDumbbells::hardCoreDiam() const 
{ return 1.0/3.0*M_PI*diameter*diameter*diameter; }


void 
IDumbbells::rescaleLengths(double scale) 
{ 
  length += scale * length;

 diameter += scale * diameter;

}

Interaction* 
IDumbbells::Clone() const 
{ return new IDumbbells(*this); }

IntEvent 
IDumbbells::getEvent(const Particle &p1,
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
      //Run this to determine when the spheres no longer intersect
      Sim->dynamics.getLiouvillean()
	.SphereSphereOutRoot(colldat, (length + diameter) * (length + diameter) ,
			     p1.testState(Particle::DYNAMIC), p2.testState(Particle::DYNAMIC));

      //colldat.dt has the upper limit of the line collision time
      //Lower limit is right now
      //Test for a line collision
      //Upper limit can be HUGE_VAL!
      if (Sim->dynamics.getLiouvillean().getOffCenterSphereOffCenterSphereCollision
	  (colldat, length, diameter, p1, p2))
	return IntEvent(p1, p2, colldat.dt, CORE, *this);
      
      return IntEvent(p1, p2, colldat.dt, WELL_OUT, *this);
    }
  else if (Sim->dynamics.getLiouvillean()
	   .SphereSphereInRoot(colldat, (length + diameter) * (length + diameter),
			       p1.testState(Particle::DYNAMIC), p2.testState(Particle::DYNAMIC))) 
    return IntEvent(p1, p2, colldat.dt, WELL_IN, *this);
  
  return IntEvent(p1, p2, HUGE_VAL, NONE, *this);
}

void
IDumbbells::runEvent(const Particle& p1, 
		  const Particle& p2,
		  const IntEvent& iEvent) const
{
  switch (iEvent.getType())
    {
    case CORE:
      {
	++Sim->eventCount;
	//We have a line interaction! Run it
	PairEventData retval(Sim->dynamics.getLiouvillean().runOffCenterSphereOffCenterSphereCollision
			     (iEvent, e, length,diameter));

	Sim->signalParticleUpdate(retval);
	
	Sim->ptrScheduler->fullUpdate(p1, p2);
	
	BOOST_FOREACH(magnet::ClonePtr<OutputPlugin> & Ptr, 
		      Sim->outputPlugins)
	  Ptr->eventUpdate(iEvent, retval);

	break;
      }
    case WELL_IN:
      {
	addToCaptureMap(p1, p2);

	//Unfortunately we cannot be smart as this well event may have
	//been pushed into both particles update lists, therefore we
	//must do a full update
	Sim->ptrScheduler->fullUpdate(p1, p2);

	Sim->freestreamAcc += iEvent.getdt();
	break;
      }
    case WELL_OUT:
      {
	removeFromCaptureMap(p1, p2);

	//Unfortunately we cannot be smart as this well event may have
	//been pushed into both particles update lists, therefore we
	//must do a full update
	Sim->ptrScheduler->fullUpdate(p1, p2);

	Sim->freestreamAcc += iEvent.getdt();
	break;
      }
    default:
      M_throw() << "Unknown collision type";
    }
}
   
void 
IDumbbells::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Type") << "Dumbbells"
      << xml::attr("Length") << length / Sim->dynamics.units().unitLength()
      << xml::attr("Elasticity") << e
      << xml::attr("Diameter") <<  diameter / Sim->dynamics.units().unitLength()
      << xml::attr("Name") << intName
      << range;

  ISingleCapture::outputCaptureMap(XML);
}

bool 
IDumbbells::captureTest(const Particle& p1, const Particle& p2) const
{
  Vector  rij = p1.getPosition() - p2.getPosition();
  Sim->dynamics.BCs().applyBC(rij);
  

  return (rij | rij) <= length*length + diameter*diameter;

}

void
IDumbbells::checkOverlaps(const Particle& part1, const Particle& part2) const
{}

void 
IDumbbells::write_povray_desc(const DYNAMO::RGB& rgb, const size_t& specID, 
			   std::ostream& os) const
{ // \todo{write this with dumbbells instead of cylinders}
  try {
    dynamic_cast<const LNOrientation&>(Sim->dynamics.getLiouvillean());
  }
  catch(std::bad_cast)
    {
      M_throw() << "Liouvillean is not an orientation liouvillean!";
    }
  
  BOOST_FOREACH(const size_t& pid, *(Sim->dynamics.getSpecies()[specID]->getRange()))
    {
      const Particle& part(Sim->particleList[pid]);

      const LNOrientation::rotData& 
	rdat(static_cast<const LNOrientation&>
	     (Sim->dynamics.getLiouvillean()).getRotData(part));

      Vector  pos(part.getPosition());
      Sim->dynamics.BCs().applyBC(pos);

      Vector  point(pos - 0.5 * length * rdat.orientation);
      
      os << "cylinder {\n <" << point[0];
      for (size_t iDim(1); iDim < NDIM; ++iDim)
	os << "," << point[iDim];

      point = pos + 0.5 * length * rdat.orientation;

      os << ">, \n <" << point[0];
      for (size_t iDim(1); iDim < NDIM; ++iDim)
	os << "," << point[iDim];

      os << ">, " << length *0.01 
	 << "\n texture { pigment { color rgb<" << rgb.R << "," << rgb.G 
	 << "," << rgb.B << "> }}\nfinish { phong 0.9 phong_size 60 }\n}\n";
    }

}
