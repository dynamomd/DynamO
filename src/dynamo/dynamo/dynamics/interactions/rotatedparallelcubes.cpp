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

#include <dynamo/dynamics/interactions/rotatedparallelcubes.hpp>
#include <dynamo/dynamics/interactions/intEvent.hpp>
#include <dynamo/dynamics/liouvillean/liouvillean.hpp>
#include <dynamo/dynamics/units/units.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <dynamo/dynamics/2particleEventData.hpp>
#include <dynamo/dynamics/BC/BC.hpp>
#include <dynamo/dynamics/ranges/1range.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/dynamics/NparticleEventData.hpp>
#include <boost/lexical_cast.hpp>
#include <magnet/xmlwriter.hpp>
#include <sstream>
#include <cmath>
#include <iomanip>

namespace dynamo {
  IParallelCubes::IParallelCubes(const magnet::xml::Node& XML, 
					       dynamo::SimData* tmp):
    Interaction(tmp,NULL)
  { operator<<(XML); }

  void 
  IParallelCubes::initialise(size_t nID)
  { 
    ID=nID; 
  }

  void 
  IParallelCubes::operator<<(const magnet::xml::Node& XML)
  { 
    if (strcmp(XML.getAttribute("Type"),"ParallelCubes"))
      M_throw() << "Attempting to load ParallelCubes from " 
		<< XML.getAttribute("Type") << " entry";
  
    Interaction::operator<<(XML);
  
    try 
      {
	_diameter = Sim->_properties.getProperty(XML.getAttribute("Diameter"),
						 Property::Units::Length());
	_e = Sim->_properties.getProperty(XML.getAttribute("Elasticity"),
					  Property::Units::Dimensionless());
	intName = XML.getAttribute("Name");
      }
    catch (boost::bad_lexical_cast &)
      {
	M_throw() << "Failed a lexical cast in CIParallelCubes";
      }
  }

  double 
  IParallelCubes::maxIntDist() const 
  { return std::sqrt(double(NDIM)) * _diameter->getMaxValue(); }

  double 
  IParallelCubes::getExcludedVolume(size_t ID) const 
  {
    double diam = _diameter->getProperty(ID);
    return diam * diam * diam; 
  }

  IntEvent 
  IParallelCubes::getEvent(const Particle &p1, const Particle &p2) const 
  { 
#ifdef DYNAMO_DEBUG
    if (!Sim->dynamics.getLiouvillean().isUpToDate(p1))
      M_throw() << "Particle 1 is not up to date";
  
    if (!Sim->dynamics.getLiouvillean().isUpToDate(p2))
      M_throw() << "Particle 2 is not up to date";
#endif

#ifdef DYNAMO_DEBUG
    if (p1 == p2)
      M_throw() << "You shouldn't pass p1==p2 events to the interactions!";
#endif 

    double d = (_diameter->getProperty(p1.getID())
		+ _diameter->getProperty(p2.getID())) * 0.5;

    double dt = Sim->dynamics.getLiouvillean().CubeCubeInRoot(p1, p2, d);

    if (dt != HUGE_VAL)
      {
#ifdef DYNAMO_OverlapTesting
	if (Sim->dynamics.getLiouvillean().cubeOverlap(p1, p2, d))
	  M_throw() << "Overlapping particles found" 
		    << ", particle1 " << p1.getID() << ", particle2 " 
		    << p2.getID() << "\nOverlap = " 
		    << Sim->dynamics.getLiouvillean().cubeOverlap(p1, p2, d)
	    / Sim->dynamics.units().unitLength();
#endif

	return IntEvent(p1, p2, dt, CORE, *this);
      }
  
    return IntEvent(p1, p2, HUGE_VAL, NONE, *this);
  }

  void
  IParallelCubes::runEvent(const Particle& p1,
				  const Particle& p2,
				  const IntEvent& iEvent) const
  {
    ++Sim->eventCount;
 
    double e = (_e->getProperty(p1.getID())
		+ _e->getProperty(p2.getID())) * 0.5;
    double d = (_diameter->getProperty(p1.getID())
		+ _diameter->getProperty(p2.getID())) * 0.5;
   
    //Run the collision and catch the data
    PairEventData EDat
      (Sim->dynamics.getLiouvillean().parallelCubeColl(iEvent, e, d)); 

    Sim->signalParticleUpdate(EDat);

    //Now we're past the event, update the scheduler and plugins
    Sim->ptrScheduler->fullUpdate(p1, p2);
  
    BOOST_FOREACH(shared_ptr<OutputPlugin> & Ptr, Sim->outputPlugins)
      Ptr->eventUpdate(iEvent,EDat);
  }
   
  void 
  IParallelCubes::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "ParallelCubes"
	<< magnet::xml::attr("Diameter") << _diameter->getName()
	<< magnet::xml::attr("Elasticity") << _e->getName()
	<< magnet::xml::attr("Name") << intName
	<< *range;
  }

  void
  IParallelCubes::checkOverlaps(const Particle& part1, const Particle& part2) const
  {
    Vector  rij = part1.getPosition() - part2.getPosition();  
    Sim->dynamics.BCs().applyBC(rij); 

    double d = (_diameter->getProperty(part1.getID())
		+ _diameter->getProperty(part2.getID())) * 0.5;
  
    if ((rij | rij) < d * d)
      derr << "Possible overlap occured in diagnostics\n ID1=" << part1.getID() 
	   << ", ID2=" << part2.getID() << "\nR_ij^2=" 
	   << (rij | rij) / pow(Sim->dynamics.units().unitLength(), 2)
	   << "\nd^2=" 
	   << d * d / pow(Sim->dynamics.units().unitLength(), 2) << std::endl;
  }
}

