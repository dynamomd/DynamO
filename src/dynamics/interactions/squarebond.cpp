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

#include "squarebond.hpp"
#include "../BC/BC.hpp"
#include "../dynamics.hpp"
#include "../units/units.hpp"
#include "../globals/global.hpp"
#include "../../simulation/particle.hpp"
#include "../interactions/intEvent.hpp"
#include "../species/species.hpp"
#include "../2particleEventData.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../../base/is_simdata.hpp"
#include "../../schedulers/scheduler.hpp"
#include "../NparticleEventData.hpp"
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <cmath>
#include <iomanip>

  
ISquareBond::ISquareBond(const magnet::xml::Node& XML, DYNAMO::SimData* tmp):
  Interaction(tmp,NULL) //A temporary value!
{ operator<<(XML); }
	    
void 
ISquareBond::operator<<(const magnet::xml::Node& XML)
{
  if (strcmp(XML.getAttribute("Type"),"SquareBond"))
    M_throw() << "Attempting to load SquareBond from non SquareBond entry";
  
  range.set_ptr(C2Range::getClass(XML,Sim));
  
  try {
    _diameter = Sim->_properties.getProperty(XML.getAttribute("Diameter"),
					     Property::Units::Length());
    _lambda = Sim->_properties.getProperty(XML.getAttribute("Lambda"),
					   Property::Units::Dimensionless());

    if (XML.getAttribute("Elasticity").valid())
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

Interaction* 
ISquareBond::Clone() const 
{ return new ISquareBond(*this); }

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
  Vector  rij = p1.getPosition() - p2.getPosition();
  Sim->dynamics.BCs().applyBC(rij);

  double d = (_diameter->getProperty(p1.getID())
	      + _diameter->getProperty(p2.getID())) * 0.5;

  double l = (_lambda->getProperty(p1.getID())
	       + _lambda->getProperty(p2.getID())) * 0.5;
  
  double ld2 = d * l * d * l;
  
#ifdef DYNAMO_DEBUG
  double d2 = d * d;
  if ((rij | rij) < d2)
    I_cerr() << "Warning! Two particles might be overlapping"
	     << "\nrij^2 = " << (rij | rij)
	     << "\nd^2 = " << d2;
#endif
  
  return (rij | rij) <= ld2;
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
    I_cerr() << std::setprecision(std::numeric_limits<float>::digits10)
	     << "Possible bonded overlap occured in diagnostics\n ID1=" << part1.getID() 
	     << ", ID2=" << part2.getID() << "\nR_ij^2=" 
	     << r2 / pow(Sim->dynamics.units().unitLength(),2)
	     << "\nd^2=" 
	     << d2 / pow(Sim->dynamics.units().unitLength(),2);
  
  if (r2 > ld2)
    I_cerr() << std::setprecision(std::numeric_limits<float>::digits10)
	     << "Possible escaped bonded pair in diagnostics\n ID1=" << part1.getID() 
	     << ", ID2=" << part2.getID() << "\nR_ij^2=" 
	     << r2 / pow(Sim->dynamics.units().unitLength(),2)
	     << "\n(lambda * d)^2=" 
	     << ld2 / pow(Sim->dynamics.units().unitLength(),2);
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

  CPDData colldat(*Sim, p1, p2);

  double d = (_diameter->getProperty(p1.getID())
	      + _diameter->getProperty(p2.getID())) * 0.5;
  double d2 = d * d;
  double l = (_lambda->getProperty(p1.getID())
	       + _lambda->getProperty(p2.getID())) * 0.5;
  
  double ld2 = d * l * d * l;


  if (Sim->dynamics.getLiouvillean()
      .SphereSphereInRoot(colldat, d2,
			  p1.testState(Particle::DYNAMIC), p2.testState(Particle::DYNAMIC)))
    {
#ifdef DYNAMO_OverlapTesting
      if (Sim->dynamics.getLiouvillean().sphereOverlap(colldat,d2))
	M_throw() << "Overlapping particles found" 
		  << ", particle1 " << p1.getID() 
		  << ", particle2 " 
		  << p2.getID() << "\nOverlap = " << (sqrt(colldat.r2) - sqrt(d2))/Sim->dynamics.units().unitLength();
#endif      
      return IntEvent(p1, p2, colldat.dt, CORE, *this);
    }
  else
    if (Sim->dynamics.getLiouvillean()
	.SphereSphereOutRoot(colldat, ld2,
			     p1.testState(Particle::DYNAMIC), p2.testState(Particle::DYNAMIC)))
      {
	return IntEvent(p1, p2, colldat.dt, BOUNCE, *this); 
      }
  
  return IntEvent(p1, p2, HUGE_VAL, NONE, *this);
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
  
  BOOST_FOREACH(magnet::ClonePtr<OutputPlugin> & Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(iEvent,EDat);

}
    
void 
ISquareBond::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Type") << "SquareBond"
      << xml::attr("Diameter") << _diameter->getName()
      << xml::attr("Lambda") << _lambda->getName()
      << xml::attr("Name") << intName
      << xml::attr("Elasticity") << _e->getName()
      << range;
}
