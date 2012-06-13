/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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

#include <dynamo/interactions/lines.hpp>
#include <dynamo/interactions/intEvent.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/2particleEventData.hpp>
#include <dynamo/BC/BC.hpp>
#include <dynamo/ranges/1range.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <cmath>
#include <iomanip>

namespace dynamo {
  ILines::ILines(const magnet::xml::Node& XML, dynamo::Simulation* tmp):
    ISingleCapture(tmp, NULL)
  {
    operator<<(XML);
  }

  void 
  ILines::initialise(size_t nID)
  {
    ID = nID; 
    ISingleCapture::initCaptureMap();
  }

  Vector ILines::getGlyphSize(size_t ID, size_t subID) const
  {
    double l = _length->getProperty(ID);
    return Vector(l, l, l);
  }

  Vector ILines::getGlyphPosition(size_t ID, size_t subID) const
  {
    Vector retval = Sim->particleList[ID].getPosition();
    Sim->BCs->applyBC(retval);
    return retval;
  }


  void 
  ILines::operator<<(const magnet::xml::Node& XML)
  { 
    if (strcmp(XML.getAttribute("Type"),"Lines"))
      M_throw() << "Attempting to load Lines from non Lines entry";
  
    Interaction::operator<<(XML);
  
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

  IntEvent 
  ILines::getEvent(const Particle &p1, const Particle &p2) const
  {
#ifdef DYNAMO_DEBUG
    if (!Sim->dynamics->isUpToDate(p1))
      M_throw() << "Particle 1 is not up to date";
  
    if (!Sim->dynamics->isUpToDate(p2))
      M_throw() << "Particle 2 is not up to date";

    if (p1 == p2)
      M_throw() << "You shouldn't pass p1==p2 events to the interactions!";
#endif 
  
    double l = (_length->getProperty(p1.getID())
		+ _length->getProperty(p2.getID())) * 0.5;

    if (isCaptured(p1, p2))
      {
	//Run this to determine when the spheres no longer intersect
	double dt = Sim->dynamics->SphereSphereOutRoot(p1, p2, l);
      
	std::pair<bool, double> colltime = Sim->dynamics->getLineLineCollision(l, p1, p2, dt);

	if (colltime.second == HUGE_VAL)
	  return IntEvent(p1, p2, dt, WELL_OUT, *this);

	//Something happens in the time interval

	if (colltime.first)
	  //Its a collision!
	  return IntEvent(p1, p2, colltime.second, CORE, *this);
	else
	  //Its a virtual event, we need to recalculate in a bit
	  return IntEvent(p1, p2, colltime.second, VIRTUAL, *this);
      }
    else 
      {
	double dt = Sim->dynamics->SphereSphereInRoot(p1, p2, l);
	if (dt != HUGE_VAL)
	  return IntEvent(p1, p2, dt, WELL_IN, *this);
      }

    return IntEvent(p1, p2, HUGE_VAL, NONE, *this);
  }

  void
  ILines::runEvent(Particle& p1, Particle& p2, const IntEvent& iEvent) const
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
	  PairEventData retval(Sim->dynamics->runLineLineCollision
			       (iEvent, e, l));

	  Sim->signalParticleUpdate(retval);
	
	  Sim->ptrScheduler->fullUpdate(p1, p2);
	
	  BOOST_FOREACH(shared_ptr<OutputPlugin> & Ptr, 
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
	<< *range;

    ISingleCapture::outputCaptureMap(XML);
  }

  bool 
  ILines::captureTest(const Particle& p1, const Particle& p2) const
  {
    if (&(*(Sim->getInteraction(p1, p2))) != this) return false;

    double l = (_length->getProperty(p1.getID())
		+ _length->getProperty(p2.getID())) * 0.5;

    return Sim->dynamics->sphereOverlap(p1, p2, l);
  }
}

