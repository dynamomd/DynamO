/*  DYNAMO:- Event driven molecular dynamics simulator 
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

#include "rotatedparallelcubes.hpp"
#include "../../dynamics/interactions/intEvent.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../units/units.hpp"
#include "../../base/is_simdata.hpp"
#include "../2particleEventData.hpp"
#include "../BC/BC.hpp"
#include "../ranges/1range.hpp"
#include "../../schedulers/scheduler.hpp"
#include "../NparticleEventData.hpp"
#include "../../datatypes/vector.xml.hpp"
#include <boost/lexical_cast.hpp>
#include <magnet/xmlwriter.hpp>
#include <sstream>
#include <cmath>
#include <iomanip>

IRotatedParallelCubes::IRotatedParallelCubes(const magnet::xml::Node& XML, 
					     DYNAMO::SimData* tmp):
  Interaction(tmp,NULL)
{ operator<<(XML); }

void 
IRotatedParallelCubes::initialise(size_t nID)
{ 
  ID=nID; 
}

void 
IRotatedParallelCubes::operator<<(const magnet::xml::Node& XML)
{ 
  if (strcmp(XML.getAttribute("Type"),"RotatedParallelCubes"))
    M_throw() << "Attempting to load RotatedParallelCubes from " 
	      << XML.getAttribute("Type") << " entry";
  
  range.set_ptr(C2Range::getClass(XML,Sim));
  
  try 
    {
      _diameter = Sim->_properties.getProperty(XML.getAttribute("Diameter"),
					       Property::Units::Length());
      _e = Sim->_properties.getProperty(XML.getAttribute("Elasticity"),
					Property::Units::Dimensionless());
      intName = XML.getAttribute("Name");
      ::operator<<(Rotation, XML.getNode("Rotation"));
    }
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in CIRotatedParallelCubes";
    }
}

double 
IRotatedParallelCubes::maxIntDist() const 
{ return std::sqrt(NDIM) * _diameter->getMaxValue(); }

double 
IRotatedParallelCubes::hardCoreDiam() const 
{ return _diameter->getMaxValue(); }


Interaction* 
IRotatedParallelCubes::Clone() const 
{ return new IRotatedParallelCubes(*this); }
  
IntEvent 
IRotatedParallelCubes::getEvent(const Particle &p1, const Particle &p2) const 
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

  CPDData colldat(*Sim, p1, p2);

  colldat.rij = Rotation * Vector(colldat.rij);
  colldat.vij = Rotation * Vector(colldat.vij);
  
  double d = (_diameter->getProperty(p1.getID())
	      + _diameter->getProperty(p2.getID())) * 0.5;

  if (Sim->dynamics.getLiouvillean().CubeCubeInRoot(colldat, d))
    {
#ifdef DYNAMO_OverlapTesting
      if (Sim->dynamics.getLiouvillean().cubeOverlap(colldat, d))
	M_throw() << "Overlapping particles found" 
		  << ", particle1 " << p1.getID() << ", particle2 " 
		  << p2.getID() << "\nOverlap = " 
		  << (sqrt(colldat.r2) - d) 
	  / Sim->dynamics.units().unitLength();
#endif

      return IntEvent(p1, p2, colldat.dt, CORE, *this);
    }
  
  return IntEvent(p1,p2,HUGE_VAL, NONE, *this);
}

void
IRotatedParallelCubes::runEvent(const Particle& p1,
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
    (Sim->dynamics.getLiouvillean().parallelCubeColl(iEvent, e, d, Rotation)); 

  Sim->signalParticleUpdate(EDat);

  //Now we're past the event, update the scheduler and plugins
  Sim->ptrScheduler->fullUpdate(p1, p2);
  
  BOOST_FOREACH(magnet::ClonePtr<OutputPlugin> & Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(iEvent,EDat);
}
   
void 
IRotatedParallelCubes::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Type") << "RotatedParallelCubes"
      << xml::attr("Diameter") << _diameter->getName()
      << xml::attr("Elasticity") << _e->getName()
      << xml::attr("Name") << intName
      << range
      << xml::tag("Rotation")
      << Rotation
      << xml::endtag("Rotation");
}

void
IRotatedParallelCubes::checkOverlaps(const Particle& part1, const Particle& part2) const
{
  Vector  rij = part1.getPosition() - part2.getPosition();  
  Sim->dynamics.BCs().applyBC(rij); 

  double d = (_diameter->getProperty(part1.getID())
	      + _diameter->getProperty(part2.getID())) * 0.5;
  
  if ((rij | rij) < d * d)
    I_cerr() << std::setprecision(std::numeric_limits<float>::digits10)
	     << "Possible overlap occured in diagnostics\n ID1=" << part1.getID() 
	     << ", ID2=" << part2.getID() << "\nR_ij^2=" 
	     << (rij | rij) / pow(Sim->dynamics.units().unitLength(), 2)
	     << "\nd^2=" 
	     << d * d / pow(Sim->dynamics.units().unitLength(), 2);
}
