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
    Interaction::initialise(nID);

    if (_et && !Sim->dynamics->hasOrientationData())
      M_throw() << "Interaction'" << getName() 
		<< "': To use a tangential coefficient of restitution, you must have orientation data for the particles in your configuration file.";
  }

  void 
  IHardSphere::operator<<(const magnet::xml::Node& XML)
  { 
    Interaction::operator<<(XML);
    _diameter = Sim->_properties.getProperty(XML.getAttribute("Diameter"), Property::Units::Length());

    if (XML.hasAttribute("Elasticity"))
      _e = Sim->_properties.getProperty(XML.getAttribute("Elasticity"), Property::Units::Dimensionless());
    
    if (XML.hasAttribute("TangentialElasticity"))
      _et = Sim->_properties.getProperty(XML.getAttribute("TangentialElasticity"), Property::Units::Dimensionless());
  }
  
  void
  IHardSphere::outputData(magnet::xml::XmlStream& XML) const
  {}

  std::array<double, 4>
  IHardSphere::getGlyphSize(size_t ID) const
  { return {{_diameter->getProperty(ID), 0, 0, 0}}; }

  double 
  IHardSphere::maxIntDist() const 
  { return _diameter->getMaxValue(); }

  double 
  IHardSphere::getExcludedVolume(size_t ID) const 
  { 
    const double diam = _diameter->getProperty(ID);
    return diam * diam * diam * M_PI / 6.0; 
  }

  Event 
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

    const double d = _diameter->getProperty(p1, p2);
    const double dt = Sim->dynamics->SphereSphereInRoot(p1, p2, d);

    if (dt != HUGE_VAL)
      return Event(p1, dt, INTERACTION, CORE, ID, p2);
  
    return Event(p1, HUGE_VAL, INTERACTION, NONE, ID, p2);
  }

  PairEventData
  IHardSphere::runEvent(Particle& p1, Particle& p2, Event iEvent)
  {
    ++Sim->eventCount;

    const double d1 = _diameter->getProperty(p1);
    const double d2 = _diameter->getProperty(p2);
    const double d = _diameter->getProperty(p1, p2);

    double e = 1.0;
    if (_e) e = _e->getProperty(p1, p2);
   
    PairEventData EDat;
    if (_et)
      return Sim->dynamics->RoughSpheresColl(iEvent, e, _et->getProperty(p1, p2), d1, d2);
    else
      return Sim->dynamics->SmoothSpheresColl(iEvent, e, d * d);
  }
   
  void 
  IHardSphere::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "HardSphere"
	<< magnet::xml::attr("Diameter") << _diameter->getName();
    if (_e) XML << magnet::xml::attr("Elasticity") << _e->getName();
    if (_et) XML << magnet::xml::attr("TangentialElasticity") << _et->getName();
    XML << magnet::xml::attr("Name") << intName
	<< range;
  }

  bool
  IHardSphere::validateState(const Particle& p1, const Particle& p2, bool textoutput) const
  {
    const double d = _diameter->getProperty(p1, p2);
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
