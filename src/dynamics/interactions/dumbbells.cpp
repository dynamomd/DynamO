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
      _length = Sim->_properties.getProperty(XML.getAttribute("Length"),
					     Property::Units::Length());
      _e = Sim->_properties.getProperty(XML.getAttribute("Elasticity"),
					       Property::Units::Dimensionless());
      _diameter = Sim->_properties.getProperty(XML.getAttribute("Diameter"),
					       Property::Units::Length());
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
{ return _length->getMaxValue() + _diameter->getMaxValue(); }

double 
IDumbbells::hardCoreDiam() const 
{ return maxIntDist(); }

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
  
  double d = (_diameter->getProperty(p1.getID())
	      + _diameter->getProperty(p2.getID())) * 0.5;

  double l = (_length->getProperty(p1.getID())
	      + _length->getProperty(p2.getID())) * 0.5;
  
  if (isCaptured(p1, p2)) 
    {
      //Run this to determine when the spheres no longer intersect
      Sim->dynamics.getLiouvillean()
	.SphereSphereOutRoot(colldat, (l + d) * (l + d),
			     p1.testState(Particle::DYNAMIC), p2.testState(Particle::DYNAMIC));

      //colldat.dt has the upper limit of the line collision time
      //Lower limit is right now
      //Test for a line collision
      //Upper limit can be HUGE_VAL!
      if (Sim->dynamics.getLiouvillean().getOffCenterSphereOffCenterSphereCollision
	  (colldat, l, d, p1, p2))
	return IntEvent(p1, p2, colldat.dt, CORE, *this);
      
      return IntEvent(p1, p2, colldat.dt, WELL_OUT, *this);
    }
  else if (Sim->dynamics.getLiouvillean()
	   .SphereSphereInRoot(colldat, (l + d) * (l + d),
			       p1.testState(Particle::DYNAMIC), 
			       p2.testState(Particle::DYNAMIC))) 
    return IntEvent(p1, p2, colldat.dt, WELL_IN, *this);
  
  return IntEvent(p1, p2, HUGE_VAL, NONE, *this);
}

void
IDumbbells::runEvent(const Particle& p1, 
		  const Particle& p2,
		  const IntEvent& iEvent) const
{
  double d = (_diameter->getProperty(p1.getID())
	      + _diameter->getProperty(p2.getID())) * 0.5;

  double l = (_length->getProperty(p1.getID())
	      + _length->getProperty(p2.getID())) * 0.5;

  double e = (_e->getProperty(p1.getID())
	      + _e->getProperty(p2.getID())) * 0.5;

  switch (iEvent.getType())
    {
    case CORE:
      {
	++Sim->eventCount;
	//We have a line interaction! Run it
	PairEventData retval(Sim->dynamics.getLiouvillean()
			     .runOffCenterSphereOffCenterSphereCollision
			     (iEvent, e, l, d));

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
      << xml::attr("Length") << _length->getName()
      << xml::attr("Elasticity") << _e->getName()
      << xml::attr("Diameter") <<  _diameter->getName()
      << xml::attr("Name") << intName
      << range;

  ISingleCapture::outputCaptureMap(XML);
}

bool 
IDumbbells::captureTest(const Particle& p1, const Particle& p2) const
{
  double d = (_diameter->getProperty(p1.getID())
	      + _diameter->getProperty(p2.getID())) * 0.5;

  double l = (_length->getProperty(p1.getID())
	      + _length->getProperty(p2.getID())) * 0.5;

  Vector  rij = p1.getPosition() - p2.getPosition();
  Sim->dynamics.BCs().applyBC(rij);

  return (rij | rij) <= (l + d) * (l + d);

}

void
IDumbbells::checkOverlaps(const Particle& part1, const Particle& part2) const
{}
