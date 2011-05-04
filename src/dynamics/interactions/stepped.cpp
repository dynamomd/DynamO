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

#include "stepped.hpp"
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
#include <boost/math/special_functions/pow.hpp>
#include <cmath>
#include <iomanip>

IStepped::IStepped(dynamo::SimData* tmp, 
		   const std::vector<steppair>& vec, C2Range* nR):
  Interaction(tmp,nR),
  _unitLength(Sim->_properties.getProperty
	      (1.0, Property::Units::Length())),
  _unitEnergy(Sim->_properties.getProperty
	      (1.0, Property::Units::Energy())),
  steps(vec)
{}

IStepped::IStepped(const magnet::xml::Node& XML, dynamo::SimData* tmp):
  Interaction(tmp, NULL) //A temporary value!
{
  operator<<(XML);
}

void 
IStepped::operator<<(const magnet::xml::Node& XML)
{
  if (strcmp(XML.getAttribute("Type"),"Stepped"))
    M_throw() << "Attempting to load Stepped from non Stepped entry";
  
  range.set_ptr(C2Range::getClass(XML,Sim));
  
  try {
    intName = XML.getAttribute("Name");

    if (!XML.getNode("Step").valid())
      M_throw() << "No steppings defined for stepped potential " 
		<< intName;

    for (magnet::xml::Node node = XML.getNode("Step"); node.valid(); ++node)
      steps.push_back(steppair(node.getAttribute("R").as<double>(),
			       node.getAttribute("E").as<double>()));
    
    std::sort(steps.rbegin(), steps.rend());

    IMultiCapture::loadCaptureMap(XML);
  }
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in CIStepped";
    }
}

Interaction* 
IStepped::Clone() const 
{ return new IStepped(*this); }

double 
IStepped::getExcludedVolume(size_t ID) const 
{ 
  //Get the inner diameter
  double diam = steps.back().first * _unitLength->getProperty(ID);
  return (M_PI / 6) * diam * diam * diam; 
}

double 
IStepped::getDiameter(size_t ID, size_t subID) const
{ return steps.back().first * _unitLength->getProperty(ID); }

Vector 
IStepped::getPosition(size_t ID, size_t subID) const
{ 
  Vector retval = Sim->particleList[ID].getPosition();
  Sim->dynamics.BCs().applyBC(retval);
  return retval;
}

double 
IStepped::maxIntDist() const 
{ return steps.front().first * _unitLength->getMaxValue(); }

void 
IStepped::initialise(size_t nID)
{
  ID = nID;
  IMultiCapture::initCaptureMap(Sim->particleList);
  
  I_cout() << "Buckets in captureMap " << captureMap.bucket_count()
	   << "\nMax bucket count " << captureMap.max_bucket_count()
	   << "\nload Factor " << captureMap.load_factor()
	   << "\nMax load Factor " << captureMap.max_load_factor();
}

int 
IStepped::captureTest(const Particle& p1, const Particle& p2) const
{
  if (&(*(Sim->dynamics.getInteraction(p1, p2))) != this) return false;
  
  Vector  rij = p1.getPosition() - p2.getPosition();
  Sim->dynamics.BCs().applyBC(rij);
  
  double r = rij.nrm();

  for (size_t i(0); i < steps.size(); ++i)
    if (r > steps[i].first * _unitLength->getMaxValue()) return i;

  return steps.size() - 1;
}

double 
IStepped::getInternalEnergy() const 
{ 
  //Once the capture maps are loaded just iterate through that determining energies
  double Energy = 0.0;

  typedef std::pair<const std::pair<size_t, size_t>, int> locpair;

  BOOST_FOREACH(const locpair& IDs, captureMap)
    Energy += steps[IDs.second - 1].second 
    * 0.5 * (_unitEnergy->getProperty(IDs.first.first)
	     + _unitEnergy->getProperty(IDs.first.second));
  
  return Energy; 
}

IntEvent
IStepped::getEvent(const Particle &p1, 
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
  const_cmap_it capstat = getCMap_it(p1,p2);

  IntEvent retval(p1, p2, HUGE_VAL, NONE, *this);

  if (capstat == captureMap.end())
    {
      double d = steps.front().first * _unitLength->getMaxValue();
      double d2 = d * d;

      //Not captured, test for capture
      if (Sim->dynamics.getLiouvillean().SphereSphereInRoot
	  (colldat, d2, p1.testState(Particle::DYNAMIC), p2.testState(Particle::DYNAMIC)))
	{
#ifdef DYNAMO_OverlapTesting
	  //Check that there is no overlap 
	  if (Sim->dynamics.getLiouvillean().sphereOverlap(colldat, d2))
	    M_throw() << "Overlapping particles found" 
		      << ", particle1 " << p1.getID() 
		      << ", particle2 " << p2.getID() 
		      << "\nOverlap = " 
		      << (sqrt(colldat.r2) - steps.front().first)
	      / Sim->dynamics.units().unitLength();
#endif
	  
	  retval = IntEvent(p1, p2, colldat.dt, WELL_IN, *this);
	}
    }
  else
    {
      //Within the potential, look for further capture or release
      //First check if there is an inner step to interact with
      if (capstat->second < static_cast<int>(steps.size()))
	{
	  double d = steps[capstat->second].first * _unitLength->getMaxValue();
	  double d2 = d * d;
	  if (Sim->dynamics.getLiouvillean().SphereSphereInRoot
	      (colldat, d2, 
	       p1.testState(Particle::DYNAMIC), p2.testState(Particle::DYNAMIC)))
	    {
#ifdef DYNAMO_OverlapTesting
	      //Check that there is no overlap 
	      if (Sim->dynamics.getLiouvillean().sphereOverlap
		  (colldat, runstepdata[capstat->second].first * l2scale))
		M_throw() << "Overlapping particles found" 
			  << ", particle1 " << p1.getID() 
			  << ", particle2 " 
			  << p2.getID() << "\nOverlap = " 
			  << (sqrt(colldat.r2) - steps[capstat->second].first)
		  /Sim->dynamics.units().unitLength();
#endif
	      
	      retval = IntEvent(p1, p2, colldat.dt, WELL_IN , *this);
	    }
	}

      {//Now test for the outward step
	double d = steps[capstat->second-1].first * _unitLength->getMaxValue();
	double d2 = d * d;
	
	if (Sim->dynamics.getLiouvillean().SphereSphereOutRoot
	    (colldat, d2, p1.testState(Particle::DYNAMIC), p2.testState(Particle::DYNAMIC)))
	  if (retval.getdt() > colldat.dt)
	    retval = IntEvent(p1, p2, colldat.dt, WELL_OUT, *this);
      }
    }

  return retval;
}

void
IStepped::runEvent(const Particle& p1, 
		    const Particle& p2,
		    const IntEvent& iEvent) const
{
  ++Sim->eventCount;

  switch (iEvent.getType())
    {
    case WELL_OUT:
      {
	cmap_it capstat = getCMap_it(p1,p2);
	
	double d = steps[capstat->second-1].first * _unitLength->getMaxValue();
	double d2 = d * d;
	double dE = steps[capstat->second-1].second;
	if (capstat->second > 1)
	  dE -= steps[capstat->second - 2].second;
	dE *= _unitEnergy->getMaxValue();

	PairEventData retVal(Sim->dynamics.getLiouvillean().SphereWellEvent
			      (iEvent, dE, d2));
	
	if (retVal.getType() != BOUNCE)
	  if (!(--capstat->second))
	    //capstat is zero so delete
	    captureMap.erase(capstat);

	Sim->signalParticleUpdate(retVal);

	Sim->ptrScheduler->fullUpdate(p1, p2);
	
	BOOST_FOREACH(magnet::ClonePtr<OutputPlugin> & Ptr, Sim->outputPlugins)
	  Ptr->eventUpdate(iEvent, retVal);
	break;
      }
    case WELL_IN:
      {
	cmap_it capstat = getCMap_it(p1, p2);
	
	if (capstat == captureMap.end())
	  capstat = captureMap.insert
	    (captureMapType::value_type
	     ((p1.getID() < p2.getID())
	      ? cMapKey(p1.getID(), p2.getID())
	      : cMapKey(p2.getID(), p1.getID()),
	      0)).first;
	
	double d = steps[capstat->second].first * _unitLength->getMaxValue();
	double d2 = d * d;
	double dE = steps[capstat->second].second;
	if (capstat->second > 0)
	  dE -= steps[capstat->second - 1].second;
	dE *= _unitEnergy->getMaxValue();


	PairEventData retVal = Sim->dynamics.getLiouvillean().SphereWellEvent(iEvent, -dE, d2);
	
	if (retVal.getType() != BOUNCE)
	  ++(capstat->second);
	else if (!capstat->second)
	  captureMap.erase(capstat);
	
	Sim->signalParticleUpdate(retVal);
	
	Sim->ptrScheduler->fullUpdate(p1, p2);
	
	BOOST_FOREACH(magnet::ClonePtr<OutputPlugin> & Ptr, Sim->outputPlugins)
	  Ptr->eventUpdate(iEvent, retVal);
	    
	break;
      }
    default:
      M_throw() << "Unknown collision type";
    } 
}

void
IStepped::checkOverlaps(const Particle& part1, const Particle& part2) const
{
  const_cmap_it capstat = getCMap_it(part1,part2);

  if (captureTest(part1,part2) != capstat->second)
    I_cerr() << "Particle " << part1.getID() << " and Particle " << part2.getID()
	     << "\nFailing as captureTest gives " << captureTest(part1,part2)
	     << "\nAnd recorded value is " << capstat->second;
}
  
void 
IStepped::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Type") << "Stepped"
      << xml::attr("Name") << intName
      << range;

  BOOST_FOREACH(const steppair& s, steps)
    XML << xml::tag("Step")
	<< xml::attr("R") << s.first 
	<< xml::attr("E") << s.second
	<< xml::endtag("Step");
  
  IMultiCapture::outputCaptureMap(XML);  
}
