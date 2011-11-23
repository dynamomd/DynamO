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

#include <dynamo/dynamics/interactions/squarebond.hpp>
#include <dynamo/dynamics/BC/BC.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/dynamics/units/units.hpp>
#include <dynamo/dynamics/globals/global.hpp>
#include <dynamo/simulation/particle.hpp>
#include <dynamo/dynamics/interactions/intEvent.hpp>
#include <dynamo/dynamics/species/species.hpp>
#include <dynamo/dynamics/2particleEventData.hpp>
#include <dynamo/dynamics/liouvillean/liouvillean.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/dynamics/NparticleEventData.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <cmath>
#include <iomanip>

namespace dynamo {
  ISquareBond::ISquareBond(const magnet::xml::Node& XML, dynamo::SimData* tmp):
    Interaction(tmp,NULL) //A temporary value!
  { operator<<(XML); }
	    
  void 
  ISquareBond::operator<<(const magnet::xml::Node& XML)
  {
    if (strcmp(XML.getAttribute("Type"),"SquareBond"))
      M_throw() << "Attempting to load SquareBond from non SquareBond entry";
  
    Interaction::operator<<(XML);
  
    try {
      _diameter = Sim->_properties.getProperty(XML.getAttribute("Diameter"),
					       Property::Units::Length());
      _lambda = Sim->_properties.getProperty(XML.getAttribute("Lambda"),
					     Property::Units::Dimensionless());

      if (XML.hasAttribute("Elasticity"))
	_e = Sim->_properties.getProperty(XML.getAttribute("Elasticity"),
					  Property::Units::Dimensionless());
      else
	_e = Sim->_properties.getProperty(1.0, Property::Units::Dimensionless());

      intName = XML.getAttribute("Name");
    }
    catch (boost::bad_lexical_cast &)
      {
	M_throw() << "Failed a lexical cast in CISquareWell";
      }
  }

  double 
  ISquareBond::getCaptureEnergy() const 
  { return 0.0; }

  double 
  ISquareBond::maxIntDist() const 
  { return _diameter->getMaxValue()
      * _lambda->getMaxValue(); }

  void 
  ISquareBond::initialise(size_t nID)
  { ID = nID; }

  bool 
  ISquareBond::captureTest(const Particle& p1, const Particle& p2) const
  {
    double d = (_diameter->getProperty(p1.getID())
		+ _diameter->getProperty(p2.getID())) * 0.5;

    double l = (_lambda->getProperty(p1.getID())
		+ _lambda->getProperty(p2.getID())) * 0.5;
  
#ifdef DYNAMO_DEBUG
    if (Sim->dynamics.getLiouvillean().sphereOverlap(p1, p2, d))
      derr << "Warning! Two particles might be overlapping"
	   << "Overlap is " << Sim->dynamics.getLiouvillean().sphereOverlap(p1, p2, d) 
	/ Sim->dynamics.units().unitLength()
	   << "\nd = " << d / Sim->dynamics.units().unitLength() << std::endl;
#endif
 
    return Sim->dynamics.getLiouvillean().sphereOverlap(p1, p2, l * d);
  }

  void
  ISquareBond::checkOverlaps(const Particle& part1, const Particle& part2) const
  {
    Vector  rij = part1.getPosition() - part2.getPosition();
    Sim->dynamics.BCs().applyBC(rij);
    double r2 = rij.nrm2();

    double d = (_diameter->getProperty(part1.getID())
		+ _diameter->getProperty(part2.getID())) * 0.5;
    double d2 = d * d;
    double l = (_lambda->getProperty(part1.getID())
		+ _lambda->getProperty(part2.getID())) * 0.5;
  
    double ld2 = d * l * d * l;


    if (r2 < d2)
      derr << "Possible bonded overlap occured in diagnostics\n ID1=" << part1.getID() 
	   << ", ID2=" << part2.getID() << "\nR_ij^2=" 
	   << r2 / pow(Sim->dynamics.units().unitLength(),2)
	   << "\nd^2=" 
	   << d2 / pow(Sim->dynamics.units().unitLength(),2) << std::endl;
  
    if (r2 > ld2)
      derr << "Possible escaped bonded pair in diagnostics\n ID1=" << part1.getID() 
	   << ", ID2=" << part2.getID() << "\nR_ij^2=" 
	   << r2 / pow(Sim->dynamics.units().unitLength(),2)
	   << "\n(lambda * d)^2=" 
	   << ld2 / pow(Sim->dynamics.units().unitLength(),2) << std::endl;
  }

  IntEvent 
  ISquareBond::getEvent(const Particle &p1, 
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

    double d = (_diameter->getProperty(p1.getID())
		+ _diameter->getProperty(p2.getID())) * 0.5;
    double l = (_lambda->getProperty(p1.getID())
		+ _lambda->getProperty(p2.getID())) * 0.5;

    IntEvent retval(p1, p2, HUGE_VAL, NONE, *this);

    double dt = Sim->dynamics.getLiouvillean()
      .SphereSphereInRoot(p1, p2, d);
    if (dt != HUGE_VAL)
      {
#ifdef DYNAMO_OverlapTesting
	if (Sim->dynamics.getLiouvillean().sphereOverlap(p1, p2, d))
	  M_throw() << "Overlapping particles found"
		    << ", particle1 " << p1.getID()
		    << ", particle2 " << p2.getID()
		    << "\nOverlap = " 
		    << Sim->dynamics.getLiouvillean()
	    .sphereOverlap(p1, p2, d)
	    / Sim->dynamics.units().unitLength();
#endif
	retval = IntEvent(p1, p2, dt, CORE, *this);
      }

    dt = Sim->dynamics.getLiouvillean().SphereSphereOutRoot(p1, p2, l * d);
    if (retval.getdt() > dt)
      retval = IntEvent(p1, p2, dt, BOUNCE, *this);
  
    return retval;
  }

  void
  ISquareBond::runEvent(const Particle& p1, const Particle& p2, 
			const IntEvent& iEvent) const
  {
    ++Sim->eventCount;

#ifdef DYNAMO_DEBUG
    if ((iEvent.getType() != BOUNCE) && (iEvent.getType() != CORE))
      M_throw() << "Unknown type found";
#endif

    double d = (_diameter->getProperty(p1.getID())
		+ _diameter->getProperty(p2.getID())) * 0.5;
    double d2 = d * d;

    double e = (_e->getProperty(p1.getID())
		+ _e->getProperty(p2.getID())) * 0.5;

    PairEventData EDat(Sim->dynamics.getLiouvillean().SmoothSpheresColl
		       (iEvent, e, d2, iEvent.getType()));

    Sim->signalParticleUpdate(EDat);
    
    //Now we're past the event, update the scheduler and plugins
    Sim->ptrScheduler->fullUpdate(p1, p2);
  
    BOOST_FOREACH(std::tr1::shared_ptr<OutputPlugin> & Ptr, Sim->outputPlugins)
      Ptr->eventUpdate(iEvent,EDat);

  }
    
  void 
  ISquareBond::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "SquareBond"
	<< magnet::xml::attr("Diameter") << _diameter->getName()
	<< magnet::xml::attr("Lambda") << _lambda->getName()
	<< magnet::xml::attr("Name") << intName
	<< magnet::xml::attr("Elasticity") << _e->getName()
	<< *range;
  }
}

