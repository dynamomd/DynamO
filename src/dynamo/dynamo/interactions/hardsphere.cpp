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

#include <dynamo/interactions/hardsphere.hpp>
#include <dynamo/interactions/intEvent.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/2particleEventData.hpp>
#include <dynamo/BC/BC.hpp>
#include <dynamo/ranges/IDRange.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/dynamics/compression.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <cmath>
#include <iomanip>

namespace dynamo {
  IHardSphere::IHardSphere(const magnet::xml::Node& XML, dynamo::Simulation* tmp):
    Interaction(tmp, NULL)
  { operator<<(XML); }

  void 
  IHardSphere::initialise(size_t nID)
  { 
    ID=nID; 
    _complete_events = 0;
    _post_event_overlap = 0;
    _accum_overlap_magnitude = 0;
    _overlapped_tests = 0;
  }

  void 
  IHardSphere::operator<<(const magnet::xml::Node& XML)
  { 
    Interaction::operator<<(XML);
    _diameter = Sim->_properties.getProperty(XML.getAttribute("Diameter"), Property::Units::Length());
    _e = Sim->_properties.getProperty(XML.getAttribute("Elasticity"), Property::Units::Dimensionless());
    intName = XML.getAttribute("Name");
  }
  
  void 
  IHardSphere::outputData(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::tag("Interaction")
	<< magnet::xml::attr("Name") << getName()
	<< magnet::xml::attr("PostEventOverlaps") << _post_event_overlap
	<< magnet::xml::attr("AvgPostEventOverlapMagnitude") << _accum_overlap_magnitude / (_post_event_overlap *  Sim->units.unitLength())
	<< magnet::xml::attr("Events") << _complete_events
	<< magnet::xml::attr("OverlapFreq") << double(_post_event_overlap) / double(_complete_events)
	<< magnet::xml::attr("OverlappedTests") << _overlapped_tests
	<< magnet::xml::endtag("Interaction");
  }

  Vector
  IHardSphere::getGlyphSize(size_t ID) const
  { 
    double diam = _diameter->getProperty(ID);
    return Vector(diam, diam, diam); 
  }

  double 
  IHardSphere::maxIntDist() const 
  { return _diameter->getMaxValue(); }

  double 
  IHardSphere::getExcludedVolume(size_t ID) const 
  { 
    double diam = _diameter->getProperty(ID);
    return diam * diam * diam * M_PI / 6.0; 
  }

  IntEvent 
  IHardSphere::getEvent(const Particle &p1, const Particle &p2) const 
  { 
#ifdef DYNAMO_DEBUG
    if (!Sim->dynamics->isUpToDate(p1))
      M_throw() << "Particle 1 is not up to date: ID1=" << p1.getID() << ", ID2=" << p2.getID() << ", delay1=" << Sim->dynamics->getParticleDelay(p1);
  
    if (!Sim->dynamics->isUpToDate(p2))
      M_throw() << "Particle 2 is not up to date: ID1=" << p1.getID() << ", ID2=" << p2.getID() << ", delay2=" << Sim->dynamics->getParticleDelay(p2);

    if (p1 == p2)
      M_throw() << "You shouldn't pass p1==p2 events to the interactions!";
#endif 

    double d = (_diameter->getProperty(p1.getID())
		 + _diameter->getProperty(p2.getID())) * 0.5;

    double dt = Sim->dynamics->SphereSphereInRoot(p1, p2, d);

    if (Sim->dynamics->sphereOverlap(p1, p2, d)) ++_overlapped_tests;

    if (dt != HUGE_VAL)
      return IntEvent(p1, p2, dt, CORE, *this);
  
    return IntEvent(p1,p2,HUGE_VAL, NONE, *this);  
  }

  void
  IHardSphere::runEvent(Particle& p1, Particle& p2, const IntEvent& iEvent)
  {
    ++Sim->eventCount;

    double d2 = (_diameter->getProperty(p1.getID())
		 + _diameter->getProperty(p2.getID())) * 0.5;
    d2 *= d2;

    const double e = (_e->getProperty(p1.getID()) + _e->getProperty(p2.getID())) * 0.5;

    PairEventData EDat(Sim->dynamics->SmoothSpheresColl(iEvent, e, d2)); 
    (*Sim->_sigParticleUpdate)(EDat);

    {
      const double d = (_diameter->getProperty(p1.getID())
		  + _diameter->getProperty(p2.getID())) * 0.5;
      const double overlap = Sim->dynamics->sphereOverlap(p1, p2, d);
      if (overlap)
	{
	  ++_post_event_overlap;
	  _accum_overlap_magnitude += overlap;
    }
}
    ++_complete_events;
    
    //Now we're past the event, update the scheduler and plugins
    Sim->ptrScheduler->fullUpdate(p1, p2);
  
    for (shared_ptr<OutputPlugin> & Ptr : Sim->outputPlugins)
      Ptr->eventUpdate(iEvent,EDat);
  }
   
  void 
  IHardSphere::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "HardSphere"
	<< magnet::xml::attr("Diameter") << _diameter->getName()
	<< magnet::xml::attr("Elasticity") << _e->getName()
	<< magnet::xml::attr("Name") << intName
	<< range;
  }

  bool
  IHardSphere::validateState(const Particle& p1, const Particle& p2, bool textoutput) const
  {
    double d = (_diameter->getProperty(p1.getID())
		 + _diameter->getProperty(p2.getID())) * 0.5;
    
    if (Sim->dynamics->sphereOverlap(p1, p2, d))
      {
	if (textoutput)
	  derr << "Particle " << p1.getID() << " and Particle " << p2.getID() 
	       << " have entered the core at " << d / Sim->units.unitLength()
	       << " and are at a distance of " 
	       << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
	       << std::endl;
	return true;
      }

    return false;
  }
}
