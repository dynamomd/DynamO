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
  IStepped::IStepped(dynamo::Simulation* tmp, shared_ptr<Potential> potential, IDPairRange* nR,
		     std::string name):
    IMultiCapture(tmp,nR),
    _unitLength(Sim->_properties.getProperty
		(Sim->units.unitLength(), 
		 Property::Units::Length())),
    _unitEnergy(Sim->_properties.getProperty
		(Sim->units.unitEnergy(), 
		 Property::Units::Energy())),
    _potential(potential)
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
  
    intName = XML.getAttribute("Name");
    
    _potential = Potential::getClass(XML.getNode("Potential"));

    IMultiCapture::loadCaptureMap(XML);
  }

  double 
  IStepped::getExcludedVolume(size_t ID) const 
  { 
    //Get the inner diameter
    double diam = _potential->hard_core_diameter() * _unitLength->getProperty(ID);
    return (M_PI / 6) * diam * diam * diam; 
  }

  Vector
  IStepped::getGlyphSize(size_t ID, size_t subID) const
  { 
    double diam = _potential->render_diameter() * _unitLength->getProperty(ID);
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
  { return (*_potential)[0].first * _unitLength->getMaxValue(); }

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

  size_t 
  IStepped::captureTest(const Particle& p1, const Particle& p2) const
  {
    if (&(*(Sim->getInteraction(p1, p2))) != this) return false;
  
    const double length_scale = 0.5 * (_unitLength->getProperty(p1.getID()) + _unitLength->getProperty(p2.getID()));

    Vector rij = p1.getPosition() - p2.getPosition();
    Sim->BCs->applyBC(rij);
    
    return _potential->calculateStepID(rij.nrm() / length_scale);
  }

  double 
  IStepped::getInternalEnergy() const 
  { 
    //Once the capture maps are loaded just iterate through that determining energies
    double Energy = 0.0;

    typedef std::pair<const std::pair<size_t, size_t>, size_t> locpair;

    BOOST_FOREACH(const locpair& IDs, captureMap)
      Energy += (*_potential)[IDs.second - 1].second
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

    const double energy_scale = 0.5 * (_unitEnergy->getProperty(p1.getID()) + _unitEnergy->getProperty(p2.getID()));
    
    return (*_potential)[capstat->second - 1].second * energy_scale;
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

    const double length_scale = 0.5 * (_unitLength->getProperty(p1.getID()) + _unitLength->getProperty(p2.getID()));

    if (capstat == captureMap.end())
      {
	double d = (*_potential)[0].first * length_scale;
	double dt = Sim->dynamics->SphereSphereInRoot(p1, p2, d);

	//Not captured, test for capture
	if (dt != HUGE_VAL)
	  retval = IntEvent(p1, p2, dt, STEP_IN, *this);
      }
    else
      {
	//Within the potential, look for further capture or release
	//First check if there is an inner step to interact with
	if (capstat->second < (*_potential).steps())
	  {
	    double d = (*_potential)[capstat->second].first * length_scale;
	    double dt = Sim->dynamics->SphereSphereInRoot(p1, p2, d);
	    
	    if (dt != HUGE_VAL)
	      retval = IntEvent(p1, p2, dt, STEP_IN , *this);
	  }

	{//Now test for the outward step
	  double d = (*_potential)[capstat->second - 1].first * length_scale;
	  double dt = Sim->dynamics->SphereSphereOutRoot(p1, p2, d);
	  if (retval.getdt() > dt)
	      retval = IntEvent(p1, p2, dt, STEP_OUT, *this);
	}
      }

    return retval;
  }

  void
  IStepped::runEvent(Particle& p1, Particle& p2, const IntEvent& iEvent) const
  {
    ++Sim->eventCount;

    const double length_scale = 0.5 * (_unitLength->getProperty(p1.getID()) + _unitLength->getProperty(p2.getID()));
    const double energy_scale = 0.5 * (_unitEnergy->getProperty(p1.getID()) + _unitEnergy->getProperty(p2.getID()));

    switch (iEvent.getType())
      {
      case STEP_OUT:
	{
	  cmap_it capstat = getCMap_it(p1,p2);
	
	  double d = (*_potential)[capstat->second-1].first * length_scale;
	  double d2 = d * d;
	  double dE = (*_potential)[capstat->second-1].second;
	  if (capstat->second > 1)
	    dE -= (*_potential)[capstat->second - 2].second;
	  dE *= energy_scale;

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
      case STEP_IN:
	{
	  cmap_it capstat = getCMap_it(p1, p2);
	
	  if (capstat == captureMap.end())
	    capstat = captureMap.insert
	      (captureMapType::value_type
	       ((p1.getID() < p2.getID())
		? cMapKey(p1.getID(), p2.getID())
		: cMapKey(p2.getID(), p1.getID()),
		0)).first;
	
	  double d = (*_potential)[capstat->second].first * length_scale;
	  double d2 = d * d;
	  double dE = (*_potential)[capstat->second].second;
	  if (capstat->second > 0)
	    dE -= (*_potential)[capstat->second - 1].second;
	  dE *= energy_scale;

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
    const size_t calculated_step_ID = captureTest(p1, p2);
    const size_t stored_step_ID = (capstat == captureMap.end()) ? 0 : capstat->second;
    const std::pair<double, double> stored_step_bounds = _potential->getStepBounds(stored_step_ID);
    const std::pair<double, double> calculated_step_bounds = _potential->getStepBounds(calculated_step_ID);
    
    if (calculated_step_ID != stored_step_ID)
      {
	if (textoutput)
	  {
	    derr << "Particle " << p1.getID() << " and Particle " << p2.getID() 
		 << " registered as being inside step " << stored_step_ID
		 << " which has limits of [" << stored_step_bounds.first
		 << ", " << stored_step_bounds.second << "] but they are at a distance of " 
		 << Sim->BCs->getDistance(p1, p2) / Sim->units.unitLength()
		 << " and this corresponds to step " << calculated_step_ID
		 << " with bounds [" << calculated_step_bounds.first << ","
		 << calculated_step_bounds.second << "]" << std::endl;
	  }
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
  
    XML << _potential;

    IMultiCapture::outputCaptureMap(XML);  
  }

  void 
  IStepped::outputData(magnet::xml::XmlStream& XML) const
  {
    using namespace magnet::xml;
    XML << tag("Interaction")
	<< attr("Name") << intName
	<< attr("Type") << "Stepped"
	<< tag("AccessedSteps");
    
    for (size_t i(0); i < _potential->cached_steps(); ++i)
      XML << tag("Step") 
	  << attr("R") << (*_potential)[i].first
	  << attr("U") << (*_potential)[i].second
	  << endtag("Step");
    
    XML << endtag("AccessedSteps")
	<< endtag("Interaction");
    
  }
}
