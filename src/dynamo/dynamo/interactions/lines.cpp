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
    Interaction::initialise(nID);
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
    _length = Sim->_properties.getProperty(XML.getAttribute("Length"), Property::Units::Length());
    _e = Sim->_properties.getProperty(XML.getAttribute("Elasticity"), Property::Units::Dimensionless());
    ICapture::loadCaptureMap(XML);   
  }

  double ILines::maxIntDist() const { return _length->getMaxValue(); }

  Event
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
  
    const double l = _length->getProperty(p1, p2);
    if (isCaptured(p1, p2))
      {
	//Run this to determine when the spheres no longer intersect
	double dt = Sim->dynamics->SphereSphereOutRoot(p1, p2, l);
      
	std::pair<bool, double> colltime = Sim->dynamics->getLineLineCollision(l, p1, p2, dt);

	if (colltime.second == HUGE_VAL)
	  return Event(p1, dt, INTERACTION, NBHOOD_OUT, ID, p2);

	//Something happens in the time interval

	if (colltime.first)
	  //Its a collision!
	  return Event(p1, colltime.second, INTERACTION, CORE, ID, p2);
	else
	  //Its a virtual event, we need to recalculate in a bit
	  return Event(p1, colltime.second, INTERACTION, VIRTUAL, ID, p2);
      }
    else 
      {
	double dt = Sim->dynamics->SphereSphereInRoot(p1, p2, l);
	if (dt != HUGE_VAL)
	  return Event(p1, dt, INTERACTION, NBHOOD_IN, ID, p2);
      }

    return Event(p1, HUGE_VAL, INTERACTION, NONE, ID, p2);
  }

  PairEventData
  ILines::runEvent(Particle& p1, Particle& p2, Event iEvent)
  {
    PairEventData retval;

    switch (iEvent._type)
      {
      case CORE:
	{
	  ++Sim->eventCount;
	  //We have a line interaction! Run it
	  const double e = _e->getProperty(p1, p2);
	  const double l = _length->getProperty(p1, p2);
	  retval = Sim->dynamics->runLineLineCollision(iEvent, e, l);
	  break;
	}
      case NBHOOD_IN:
	{
	  ICapture::add(p1, p2);
	  retval = PairEventData(p1, p2, *Sim->species(p1), *Sim->species(p2), VIRTUAL);
	  iEvent._type = VIRTUAL;
	  break;
	}
      case NBHOOD_OUT:
	{
	  ICapture::remove(p1, p2);
	  retval = PairEventData(p1, p2, *Sim->species(p1), *Sim->species(p2), VIRTUAL);
	  iEvent._type = VIRTUAL;
	  break;
	}
      case VIRTUAL:
	{
	  retval = PairEventData(p1, p2, *Sim->species(p1), *Sim->species(p2), VIRTUAL);
	  iEvent._type = VIRTUAL;
	  break;
	}
      default:
	M_throw() << "Unknown collision type";
      }

    return retval;
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

    const double l = _length->getProperty(p1, p2);
    return Sim->dynamics->sphereOverlap(p1, p2, l) > 0;
  }

  bool
  ILines::validateState(const Particle& p1, const Particle& p2, bool textoutput) const
  {
    const double l = _length->getProperty(p1, p2);

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
