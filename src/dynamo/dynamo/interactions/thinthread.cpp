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

#include <dynamo/interactions/thinthread.hpp>
#include <dynamo/BC/BC.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/globals/global.hpp>
#include <dynamo/particle.hpp>
#include <dynamo/interactions/intEvent.hpp>
#include <dynamo/species/species.hpp>
#include <dynamo/2particleEventData.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <cmath>
#include <iomanip>

namespace dynamo {
  IThinThread::IThinThread(const magnet::xml::Node& XML, dynamo::Simulation* tmp):
    ISquareWell(tmp, NULL) //A temporary value!
  {
    operator<<(XML);
  }

  void 
  IThinThread::operator<<(const magnet::xml::Node& XML)
  {
    Interaction::operator<<(XML);
    _diameter = Sim->_properties.getProperty(XML.getAttribute("Diameter"), Property::Units::Length());
    _lambda = Sim->_properties.getProperty(XML.getAttribute("Lambda"), Property::Units::Dimensionless());
    _wellDepth = Sim->_properties.getProperty(XML.getAttribute("WellDepth"), Property::Units::Energy());
    if (XML.hasAttribute("Elasticity"))
      _e = Sim->_properties.getProperty(XML.getAttribute("Elasticity"), Property::Units::Dimensionless());
    else
      _e = Sim->_properties.getProperty(1.0, Property::Units::Dimensionless());
    ICapture::loadCaptureMap(XML);   
  }

  IntEvent
  IThinThread::getEvent(const Particle &p1, 
			const Particle &p2) const 
  {
#ifdef DYNAMO_DEBUG
    if (!Sim->dynamics->isUpToDate(p1))
      M_throw() << "Particle 1 is not up to date";
  
    if (!Sim->dynamics->isUpToDate(p2))
      M_throw() << "Particle 2 is not up to date";

    if (p1 == p2)
      M_throw() << "You shouldn't pass p1==p2 events to the interactions!";
#endif 

    const double d = _diameter->getProperty(p1, p2);
    const double l = _lambda->getProperty(p1, p2);

    IntEvent retval(p1, p2, HUGE_VAL, NONE, *this);

    if (isCaptured(p1, p2))
      {
	double dt = Sim->dynamics->SphereSphereInRoot(p1, p2, d);
	if (dt != HUGE_VAL)
	  retval = IntEvent(p1, p2, dt, CORE, *this);

	dt = Sim->dynamics->SphereSphereOutRoot(p1, p2, l * d);
	if (retval.getdt() > dt)
	  retval = IntEvent(p1, p2, dt, STEP_OUT, *this);
      }
    else
      {
	double dt = Sim->dynamics->SphereSphereInRoot(p1, p2, d);

	if (dt != HUGE_VAL)
	  retval = IntEvent(p1, p2, dt, CORE, *this);
      }

    return retval;
  }

  PairEventData
  IThinThread::runEvent(Particle& p1, Particle& p2, const IntEvent& iEvent)
  {
    ++Sim->eventCount;

    const double d = _diameter->getProperty(p1, p2);
    const double d2 = d * d;

    const double e = _e->getProperty(p1, p2);

    const double l = _lambda->getProperty(p1, p2);
    const double ld2 = d * l * d * l;

    const double wd = _wellDepth->getProperty(p1, p2);

    PairEventData retVal;
    switch (iEvent.getType())
      {
      case CORE:
	{
	  retVal = Sim->dynamics->SmoothSpheresColl(iEvent, e, d2, CORE);
	  if (!isCaptured(p1, p2))
	    {
	      retVal.setType(STEP_IN);
	      ICapture::add(p1, p2);
	    }
	  break;
	}
      case STEP_OUT:
	{
	  retVal = Sim->dynamics->SphereWellEvent(iEvent, -wd, ld2, 0);
	  if (retVal.getType() != BOUNCE) ICapture::remove(p1, p2);
	  break;
	}
      default:
	M_throw() << "Unknown collision type";
      } 

    return retVal;
  }

  bool
  IThinThread::validateState(const Particle& p1, const Particle& p2, bool textoutput) const
  {
    const double d = _diameter->getProperty(p1, p2);
    const double l = _lambda->getProperty(p1, p2);

    if (isCaptured(p1, p2))
      {
	if (!Sim->dynamics->sphereOverlap(p1, p2, l * d))
	  {
	    if (textoutput)
	      derr << "Particle " << p1.getID() << " and Particle " << p2.getID() 
		   << " registered as being inside the well at " << l * d / Sim->units.unitLength()
		   << " but they are at a distance of " 
		   << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
		   << std::endl;
	    
	    return true;
	  }

	if (Sim->dynamics->sphereOverlap(p1, p2, d))
	  {
	    if (textoutput)
	      derr << "Particle " << p1.getID() << " and Particle " << p2.getID() 
		   << " are inside the well with an inner hard core at " << d / Sim->units.unitLength()
		   << " but they are at a distance of " 
		   << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
		   << std::endl;
	    
	    return true;
	  }
      }
    else
      if (Sim->dynamics->sphereOverlap(p1, p2, d))
	{
	  if (textoutput)
	    derr << "Particle " << p1.getID() << " and Particle " << p2.getID() 
		 << " have entered the core at " << d / Sim->units.unitLength()
		 << " and are at a distance of "
		 << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
		 << " AND they've not entered the thin-thread well either."
		 << std::endl;
	  
	  return true;
	}

    return false;
  }

  void 
  IThinThread::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "ThinThread"
	<< magnet::xml::attr("Diameter") << _diameter->getName()
	<< magnet::xml::attr("Elasticity") << _e->getName()
	<< magnet::xml::attr("Lambda") << _lambda->getName()
	<< magnet::xml::attr("WellDepth") << _wellDepth->getName()
	<< magnet::xml::attr("Name") << intName
	<< *range;
  
    ICapture::outputCaptureMap(XML);  
  }
}
