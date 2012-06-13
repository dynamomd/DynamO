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

#include <dynamo/interactions/roughhardsphere.hpp>
#include <dynamo/interactions/intEvent.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/2particleEventData.hpp>
#include <dynamo/BC/BC.hpp>
#include <dynamo/ranges/1range.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/dynamics/compression.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>
#include <cmath>
#include <iomanip>

namespace dynamo {
  IRoughHardSphere::IRoughHardSphere(const magnet::xml::Node& XML, dynamo::Simulation* tmp):
    Interaction(tmp,NULL)
  {
    operator<<(XML);
  }

  void 
  IRoughHardSphere::initialise(size_t nID)
  { ID=nID; }

  void 
  IRoughHardSphere::operator<<(const magnet::xml::Node& XML)
  { 
    if (strcmp(XML.getAttribute("Type"),"RoughHardSphere"))
      M_throw() << "Attempting to load Hardsphere from non RoughHardSphere entry";
  
    Interaction::operator<<(XML);
  
    try 
      {
	_diameter = Sim->_properties.getProperty(XML.getAttribute("Diameter"),
						 Property::Units::Length());
	_e = Sim->_properties.getProperty(XML.getAttribute("Elasticity"),
					  Property::Units::Dimensionless());
	_et = Sim->_properties.getProperty(XML.getAttribute("TangentialElasticity"),
					   Property::Units::Dimensionless());
	intName = XML.getAttribute("Name");
      }
    catch (boost::bad_lexical_cast &)
      {
	M_throw() << "Failed a lexical cast in IRoughHardSphere";
      }
  }

  double 
  IRoughHardSphere::maxIntDist() const 
  { return _diameter->getMaxValue(); }

  double 
  IRoughHardSphere::getExcludedVolume(size_t ID) const 
  { 
    double diam = _diameter->getProperty(ID);
    return diam * diam * diam * M_PI / 6.0; 
  }

  Vector 
  IRoughHardSphere::getGlyphSize(size_t ID, size_t subID) const
  { 
    double diam = _diameter->getProperty(ID);
    return Vector(diam, diam, diam); 
  }

  Vector 
  IRoughHardSphere::getGlyphPosition(size_t ID, size_t subID) const
  { 
    Vector retval = Sim->particleList[ID].getPosition();
    Sim->BCs->applyBC(retval);
    return retval;
  }

  IntEvent 
  IRoughHardSphere::getEvent(const Particle& p1, const Particle& p2) const 
  { 
#ifdef DYNAMO_DEBUG
    if (!Sim->dynamics->isUpToDate(p1))
      M_throw() << "Particle 1 is not up to date";
  
    if (!Sim->dynamics->isUpToDate(p2))
      M_throw() << "Particle 2 is not up to date";
#endif

#ifdef DYNAMO_DEBUG
    if (p1 == p2)
      M_throw() << "You shouldn't pass p1==p2 events to the interactions!";
#endif 

    double d = (_diameter->getProperty(p1.getID())
		+ _diameter->getProperty(p2.getID())) * 0.5;

    double dt = Sim->dynamics->SphereSphereInRoot(p1, p2, d);
    if (dt != HUGE_VAL)
      {
#ifdef DYNAMO_OverlapTesting
	if (Sim->dynamics->sphereOverlap(p1, p2, d))
	  M_throw() << "Overlapping particles found"
		    << ", particle1 " << p1.getID()
		    << ", particle2 " << p2.getID()
		    << "\nOverlap = " 
		    << Sim->dynamics.getDynamics()
	    .sphereOverlap(p1, p2, d)
	    / Sim->units.unitLength();
#endif

	return IntEvent(p1, p2, dt, CORE, *this);
      }
  
    return IntEvent(p1, p2, HUGE_VAL, NONE, *this);
  }

  void
  IRoughHardSphere::runEvent(Particle& p1, Particle& p2, const IntEvent& iEvent) const
  {
    ++Sim->eventCount;
    
    double e = (_e->getProperty(p1.getID())
		+ _e->getProperty(p2.getID())) * 0.5;

    double et = (_et->getProperty(p1.getID())
		 + _et->getProperty(p2.getID())) * 0.5;

    double d2 = (_diameter->getProperty(p1.getID())
		 + _diameter->getProperty(p2.getID())) * 0.5;
    d2 *= d2;

    //Run the collision and catch the data
    PairEventData EDat
      (Sim->dynamics->RoughSpheresColl(iEvent, e, et, d2)); 

    Sim->signalParticleUpdate(EDat);

    //Now we're past the event, update the scheduler and plugins
    Sim->ptrScheduler->fullUpdate(p1, p2);
  
    BOOST_FOREACH(shared_ptr<OutputPlugin> & Ptr, Sim->outputPlugins)
      Ptr->eventUpdate(iEvent,EDat);
  }
   
  void 
  IRoughHardSphere::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "HardSphere"
	<< magnet::xml::attr("Diameter") << _diameter->getName()
	<< magnet::xml::attr("Elasticity") << _e->getName()
	<< magnet::xml::attr("TangentialElasticity") << _et->getName()
	<< magnet::xml::attr("Name") << intName
	<< *range;
  }

  void
  IRoughHardSphere::checkOverlaps(const Particle& part1, const Particle& part2) const
  {
    Vector  rij = part1.getPosition() - part2.getPosition();  
    Sim->BCs->applyBC(rij); 

    double d2 = (_diameter->getProperty(part1.getID())
		 + _diameter->getProperty(part2.getID())) * 0.5;
    d2 *= d2;
  
    if ((rij | rij) < d2)
      derr << "Possible overlap occured in diagnostics\n ID1=" << part1.getID() 
	   << ", ID2=" << part2.getID() << "\nR_ij^2=" 
	   << (rij | rij) / pow(Sim->units.unitLength(),2)
	   << "\nd^2=" 
	   << d2 / pow(Sim->units.unitLength(),2) << std::endl;
  }
}

