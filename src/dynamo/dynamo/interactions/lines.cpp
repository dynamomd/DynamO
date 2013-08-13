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
#include <dynamo/ranges/IDRange.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <cmath>
#include <iomanip>

namespace dynamo {
  ILines::ILines(const magnet::xml::Node& XML, dynamo::Simulation* tmp):
    ICapture(tmp, NULL)
  {
    operator<<(XML);
  }

  void 
  ILines::initialise(size_t nID)
  {
    ID = nID; 
    ICapture::initCaptureMap();
  }

  std::array<double, 4> ILines::getGlyphSize(size_t ID) const
  {
    return {{_length->getProperty(ID), 0, 0, 0}};
  }

  void 
  ILines::operator<<(const magnet::xml::Node& XML)
  { 
    Interaction::operator<<(XML);
    _length = Sim->_properties.getProperty(XML.getAttribute("Length"),
					   Property::Units::Length());
    _e = Sim->_properties.getProperty(XML.getAttribute("Elasticity"),
				      Property::Units::Dimensionless());
    intName = XML.getAttribute("Name");
    ICapture::loadCaptureMap(XML);   
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
	  return IntEvent(p1, p2, dt, NBHOOD_OUT, *this);

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
	  return IntEvent(p1, p2, dt, NBHOOD_IN, *this);
      }

    return IntEvent(p1, p2, HUGE_VAL, NONE, *this);
  }

  void
  ILines::runEvent(Particle& p1, Particle& p2, const IntEvent& iEvent)
  {
    PairEventData retval;

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

	  retval = Sim->dynamics->runLineLineCollision(iEvent, e, l);
	  break;
	}
      case NBHOOD_IN:
	{
	  ICapture::add(p1, p2);
	  retval = PairEventData(p1, p2, *Sim->species[p1], *Sim->species[p2], VIRTUAL);
	  iEvent.setType(VIRTUAL);
	  break;
	}
      case NBHOOD_OUT:
	{
	  ICapture::remove(p1, p2);
	  retval = PairEventData(p1, p2, *Sim->species[p1], *Sim->species[p2], VIRTUAL);
	  iEvent.setType(VIRTUAL);
	  break;
	}
      case VIRTUAL:
	{
	  retval = PairEventData(p1, p2, *Sim->species[p1], *Sim->species[p2], VIRTUAL);
	  iEvent.setType(VIRTUAL);
	  break;
	}
      default:
	M_throw() << "Unknown collision type";
      }
    
    Sim->_sigParticleUpdate(retval);
    
    Sim->ptrScheduler->fullUpdate(p1, p2);
    
    for (shared_ptr<OutputPlugin> & Ptr : Sim->outputPlugins)
      Ptr->eventUpdate(iEvent, retval);
  }
   
  void 
  ILines::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "Lines"
	<< magnet::xml::attr("Length") << _length->getName()
	<< magnet::xml::attr("Elasticity") << _e->getName()
	<< magnet::xml::attr("Name") << intName
	<< range;

    ICapture::outputCaptureMap(XML);
  }

  size_t
  ILines::captureTest(const Particle& p1, const Particle& p2) const
  {
    if (&(*(Sim->getInteraction(p1, p2))) != this) return false;

    double l = (_length->getProperty(p1.getID())
		+ _length->getProperty(p2.getID())) * 0.5;

    return Sim->dynamics->sphereOverlap(p1, p2, l) > 0;
  }

  bool
  ILines::validateState(const Particle& p1, const Particle& p2, bool textoutput) const
  {
    double l = (_length->getProperty(p1.getID())
		+ _length->getProperty(p2.getID())) * 0.5;

    if (isCaptured(p1, p2))
      {
	if (!Sim->dynamics->sphereOverlap(p1, p2, l))
	  {
	    if (textoutput)
	      derr << "Particle " << p1.getID() << " and Particle " << p2.getID() 
		   << " are registered as being closer than " << l / Sim->units.unitLength()
		   << " but they are at a distance of " 
		   << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
		   << std::endl;

	    return true;
	  }
      }
    else
      if (Sim->dynamics->sphereOverlap(p1, p2, l))
	{
	  if (textoutput)
	    derr << "Particle " << p1.getID() << " and Particle " << p2.getID() 
		 << " are not registered as being closer than " << l / Sim->units.unitLength()
		 << " but they are at a distance of " 
		 << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
		 << std::endl;
	  return true;
	}

    return false;
  }
}
