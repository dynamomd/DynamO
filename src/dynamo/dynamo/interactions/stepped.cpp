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

#include <dynamo/interactions/stepped.hpp>
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
#include <boost/math/special_functions/pow.hpp>
#include <cmath>
#include <iomanip>

namespace dynamo {
  IStepped::IStepped(dynamo::Simulation* tmp, 
		     const std::vector<steppair>& vec, IDPairRange* nR,
		     std::string name):
    IMultiCapture(tmp,nR),
    _unitLength(Sim->_properties.getProperty
		(Sim->units.unitLength(), 
		 Property::Units::Length())),
    _unitEnergy(Sim->_properties.getProperty
		(Sim->units.unitEnergy(), 
		 Property::Units::Energy())),
    steps(vec)
  { intName = name; }

  IStepped::IStepped(const magnet::xml::Node& XML, dynamo::Simulation* tmp):
    IMultiCapture(tmp, NULL), //A temporary value!
    _unitLength(Sim->_properties.getProperty
		(Sim->units.unitLength(), 
		 Property::Units::Length())),
    _unitEnergy(Sim->_properties.getProperty
		(Sim->units.unitEnergy(), 
		 Property::Units::Energy()))
  {
    operator<<(XML);
  }

  void 
  IStepped::operator<<(const magnet::xml::Node& XML)
  {
    Interaction::operator<<(XML);
  
    try {
      intName = XML.getAttribute("Name");

      if (!XML.hasNode("Step"))
	M_throw() << "No steppings defined for stepped potential " 
		  << intName;

      for (magnet::xml::Node node = XML.fastGetNode("Step"); node.valid(); ++node)
	steps.push_back(steppair(node.getAttribute("R").as<double>(),
				 node.getAttribute("E").as<double>()));
    
      std::sort(steps.rbegin(), steps.rend());

      IMultiCapture::loadCaptureMap(XML);
    }
    catch (boost::bad_lexical_cast &)
      {
	M_throw() << "Failed a lexical cast in CIStepped";
      }

    if (steps.empty())
      M_throw() << "No steps defined in SteppedPotential Interaction with name " 
		<< getName();
  }

  double 
  IStepped::getExcludedVolume(size_t ID) const 
  { 
    //Get the inner diameter
    double diam = steps.back().first * _unitLength->getProperty(ID);
    return (M_PI / 6) * diam * diam * diam; 
  }

  Vector
  IStepped::getGlyphSize(size_t ID, size_t subID) const
  { 
    double diam = steps.back().first * _unitLength->getProperty(ID);
    return Vector(diam, diam, diam);
  }

  Vector 
  IStepped::getGlyphPosition(size_t ID, size_t subID) const
  { 
    Vector retval = Sim->particles[ID].getPosition();
    Sim->BCs->applyBC(retval);
    return retval;
  }

  double 
  IStepped::maxIntDist() const 
  { return steps.front().first * _unitLength->getMaxValue(); }

  void 
  IStepped::initialise(size_t nID)
  {
    ID = nID;
    IMultiCapture::initCaptureMap();
  
    dout << "Buckets in captureMap " << captureMap.bucket_count()
	 << "\nMax bucket count " << captureMap.max_bucket_count()
	 << "\nload Factor " << captureMap.load_factor()
	 << "\nMax load Factor " << captureMap.max_load_factor() << std::endl;
  }

  int 
  IStepped::captureTest(const Particle& p1, const Particle& p2) const
  {
    if (&(*(Sim->getInteraction(p1, p2))) != this) return false;
  
    Vector  rij = p1.getPosition() - p2.getPosition();
    Sim->BCs->applyBC(rij);
  
    double r = rij.nrm();

    //Uncaptured value
    size_t retval = 0;

    //Check when it is less
    for (size_t i(0); i < steps.size(); ++i)
      {
	if (r > steps[i].first * _unitLength->getMaxValue()) 
	  break;
	retval = i+1;
      }

    return retval;
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

  double 
  IStepped::getInternalEnergy(const Particle& p1, const Particle& p2) const
  {
    const_cmap_it capstat = getCMap_it(p1,p2);
    if (capstat == captureMap.end())
      return 0;
    else
      return steps[capstat->second - 1].second
	* 0.5 * (_unitEnergy->getProperty(p1.getID())
		 + _unitEnergy->getProperty(p2.getID()))
	* isCaptured(p1, p2);
  }

  IntEvent
  IStepped::getEvent(const Particle &p1, 
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

    const_cmap_it capstat = getCMap_it(p1,p2);

    IntEvent retval(p1, p2, HUGE_VAL, NONE, *this);

    if (capstat == captureMap.end())
      {
	double d = steps.front().first * _unitLength->getMaxValue();
	double dt 
	  = Sim->dynamics->SphereSphereInRoot(p1, p2, d);

	//Not captured, test for capture
	if (dt != HUGE_VAL)
	  retval = IntEvent(p1, p2, dt, WELL_IN, *this);
      }
    else
      {
	//Within the potential, look for further capture or release
	//First check if there is an inner step to interact with
	if (capstat->second < static_cast<int>(steps.size()))
	  {
	    double d = steps[capstat->second].first * _unitLength->getMaxValue();
	    double dt = Sim->dynamics->SphereSphereInRoot
	      (p1, p2, d);
	    
	    if (dt != HUGE_VAL)
	      retval = IntEvent(p1, p2, dt, WELL_IN , *this);
	  }

	{//Now test for the outward step
	  double d = steps[capstat->second-1].first * _unitLength->getMaxValue();
	  double dt = Sim->dynamics->SphereSphereOutRoot(p1, p2, d);
	  if (retval.getdt() > dt)
	      retval = IntEvent(p1, p2, dt, WELL_OUT, *this);
	}
      }

    return retval;
  }

  void
  IStepped::runEvent(Particle& p1, Particle& p2, const IntEvent& iEvent) const
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

	  PairEventData retVal(Sim->dynamics->SphereWellEvent
			       (iEvent, dE, d2));
	
	  if (retVal.getType() != BOUNCE)
	    if (!(--capstat->second))
	      //capstat is zero so delete
	      captureMap.erase(capstat);

	  Sim->signalParticleUpdate(retVal);

	  Sim->ptrScheduler->fullUpdate(p1, p2);
	
	  BOOST_FOREACH(shared_ptr<OutputPlugin> & Ptr, Sim->outputPlugins)
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


	  PairEventData retVal = Sim->dynamics->SphereWellEvent(iEvent, -dE, d2);
	
	  if (retVal.getType() != BOUNCE)
	    ++(capstat->second);
	  else if (!capstat->second)
	    captureMap.erase(capstat);
	
	  Sim->signalParticleUpdate(retVal);
	
	  Sim->ptrScheduler->fullUpdate(p1, p2);
	
	  BOOST_FOREACH(shared_ptr<OutputPlugin> & Ptr, Sim->outputPlugins)
	    Ptr->eventUpdate(iEvent, retVal);
	    
	  break;
	}
      default:
	M_throw() << "Unknown collision type";
      } 
  }

  bool
  IStepped::validateState(const Particle& p1, const Particle& p2, bool textoutput) const
  {
    const_cmap_it capstat = getCMap_it(p1, p2);
    int val = captureTest(p1, p2);

    if (capstat == captureMap.end())
      {
	if (val != 0)
	  {
	    double d = steps.front().first * _unitLength->getMaxValue();

	    if (textoutput)
	      derr << "Particle " << p1.getID() << " and Particle " << p2.getID() 
		   << " registered as being outside the steps, starting at " << d / Sim->units.unitLength()
		   << " but they are at a distance of "
		   << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
		   << std::endl;

	    return true;
	  }
      }
    else
      if (capstat->second != val)
	{
	  if (textoutput)
	    derr << "Particle " << p1.getID() << " and Particle " << p2.getID() 
		 << " registered as being inside step " << capstat->second 
		 << " which has limits of [" << ((capstat->second < static_cast<int>(steps.size())) ? 
						 steps[capstat->second].first * _unitLength->getMaxValue() : 0) 
	      / Sim->units.unitLength()
		 << ", " << steps[capstat->second-1].first * _unitLength->getMaxValue() / Sim->units.unitLength()
		 << "] but they are at a distance of " 
		 << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
		 << " and this corresponds to step " << val
		 << std::endl;
	    
	  return true;
	}

    return false;
  }

  void 
  IStepped::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "Stepped"
	<< magnet::xml::attr("Name") << intName
	<< *range;

    BOOST_FOREACH(const steppair& s, steps)
      XML << magnet::xml::tag("Step")
	  << magnet::xml::attr("R") << s.first 
	  << magnet::xml::attr("E") << s.second
	  << magnet::xml::endtag("Step");
  
    IMultiCapture::outputCaptureMap(XML);  
  }
}
