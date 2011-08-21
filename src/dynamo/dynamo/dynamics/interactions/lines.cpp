/*  dynamo:- Event driven molecular dynamics simulator 
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

#include <dynamo/dynamics/interactions/lines.hpp>
#include <dynamo/dynamics/interactions/intEvent.hpp>
#include <dynamo/dynamics/liouvillean/liouvillean.hpp>
#include <dynamo/dynamics/units/units.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <dynamo/dynamics/2particleEventData.hpp>
#include <dynamo/dynamics/BC/BC.hpp>
#include <dynamo/dynamics/ranges/1range.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/dynamics/NparticleEventData.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <cmath>
#include <iomanip>

namespace dynamo {
  ILines::ILines(const magnet::xml::Node& XML, dynamo::SimData* tmp):
    Interaction(tmp, NULL)
  {
    operator<<(XML);
  }

  void 
  ILines::initialise(size_t nID)
  {
    ID = nID; 
  
    ISingleCapture::initCaptureMap(Sim->particleList);
  }

  void 
  ILines::operator<<(const magnet::xml::Node& XML)
  { 
    if (strcmp(XML.getAttribute("Type"),"Lines"))
      M_throw() << "Attempting to load Lines from non Lines entry";
  
    range.set_ptr(C2Range::getClass(XML,Sim));
  
    try 
      {
	_length = Sim->_properties.getProperty(XML.getAttribute("Length"),
					       Property::Units::Length());
	_e = Sim->_properties.getProperty(XML.getAttribute("Elasticity"),
					  Property::Units::Dimensionless());
	intName = XML.getAttribute("Name");
	ISingleCapture::loadCaptureMap(XML);   
      }
    catch (boost::bad_lexical_cast &)
      {
	M_throw() << "Failed a lexical cast in CILines";
      }
  }

  double 
  ILines::maxIntDist() const 
  { return _length->getMaxValue(); }

  Interaction* 
  ILines::Clone() const 
  { return new ILines(*this); }

  IntEvent 
  ILines::getEvent(const Particle &p1,
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
  
    double l = (_length->getProperty(p1.getID())
		+ _length->getProperty(p2.getID())) * 0.5;
    double l2 = l*l;

    if (isCaptured(p1, p2)) 
      {
	//Run this to determine when the spheres no longer intersect
	Sim->dynamics.getLiouvillean()
	  .SphereSphereOutRoot(colldat, l2,
			       p1.testState(Particle::DYNAMIC), p2.testState(Particle::DYNAMIC));
      
	//colldat.dt has the upper limit of the line collision time
	//Lower limit is right now
	//Test for a line collision
	//Upper limit can be HUGE_VAL!
	if (Sim->dynamics.getLiouvillean().getLineLineCollision
	    (colldat, l, p1, p2))
	  return IntEvent(p1, p2, colldat.dt, CORE, *this);
      
	return IntEvent(p1, p2, colldat.dt, WELL_OUT, *this);
      }
    else if (Sim->dynamics.getLiouvillean()
	     .SphereSphereInRoot(colldat, l2,
				 p1.testState(Particle::DYNAMIC), p2.testState(Particle::DYNAMIC))) 
      return IntEvent(p1, p2, colldat.dt, WELL_IN, *this);
  
    return IntEvent(p1, p2, HUGE_VAL, NONE, *this);
  }

  void
  ILines::runEvent(const Particle& p1, 
		   const Particle& p2,
		   const IntEvent& iEvent) const
  {
    switch (iEvent.getType())
      {
      case CORE:
	{
	  ++Sim->eventCount;
	  //We have a line interaction! Run it
	  double e = (_e->getProperty(p1.getID())
		      + _e->getProperty(p2.getID())) * 0.5;
	  double l = (_length->getProperty(p1.getID())
		      + _length->getProperty(p2.getID())) * 0.5;
	  PairEventData retval(Sim->dynamics.getLiouvillean().runLineLineCollision
			       (iEvent, e, l));

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
  ILines::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "Lines"
	<< magnet::xml::attr("Length") << _length->getName()
	<< magnet::xml::attr("Elasticity") << _e->getName()
	<< magnet::xml::attr("Name") << intName
	<< range;

    ISingleCapture::outputCaptureMap(XML);
  }

  bool 
  ILines::captureTest(const Particle& p1, const Particle& p2) const
  {
    if (&(*(Sim->dynamics.getInteraction(p1, p2))) != this) return false;

    Vector  rij = p1.getPosition() - p2.getPosition();
    Sim->dynamics.BCs().applyBC(rij);
 
    double l2 = (_length->getProperty(p1.getID())
		 + _length->getProperty(p2.getID())) * 0.5;
    l2 *= l2;
 
    return (rij | rij) <= l2;
  }
}

