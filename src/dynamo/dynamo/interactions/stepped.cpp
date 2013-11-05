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
#include <dynamo/outputplugins/tickerproperty/radialdist.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <cmath>
#include <iomanip>

namespace dynamo {
  IStepped::IStepped(const magnet::xml::Node& XML, dynamo::Simulation* tmp):
    ICapture(tmp, NULL),
    _lengthScale(Sim->_properties.getProperty(Sim->units.unitLength(), Property::Units::Length())),
    _energyScale(Sim->_properties.getProperty(Sim->units.unitEnergy(), Property::Units::Energy()))
  {
    operator<<(XML);
  }

  void 
  IStepped::operator<<(const magnet::xml::Node& XML)
  {
    Interaction::operator<<(XML);
  
    intName = XML.getAttribute("Name");
    
    _potential = Potential::getClass(XML.getNode("Potential"));

    _lengthScale = Sim->_properties.getProperty(XML.getAttribute("LengthScale"), Property::Units::Length());

    _energyScale = Sim->_properties.getProperty(XML.getAttribute("EnergyScale"), Property::Units::Energy());

    ICapture::loadCaptureMap(XML);
  }

  double 
  IStepped::getExcludedVolume(size_t ID) const 
  { 
    //Get the inner diameter
    double diam = _potential->hard_core_diameter() * _lengthScale->getProperty(ID);
    return (M_PI / 6) * diam * diam * diam; 
  }

  std::array<double, 4>
  IStepped::getGlyphSize(size_t ID) const
  { 
    return {{_potential->render_diameter() * _lengthScale->getProperty(ID), 0, 0, 0}};
  }

  double 
  IStepped::maxIntDist() const 
  { return _potential->max_distance() * _lengthScale->getMaxValue(); }

  void 
  IStepped::initialise(size_t nID)
  {
    Interaction::initialise(nID);
    ICapture::initCaptureMap();
  }

  size_t 
  IStepped::captureTest(const Particle& p1, const Particle& p2) const
  {
    if (&(*(Sim->getInteraction(p1, p2))) != this) return false;
  
    const double length_scale = _lengthScale->getProperty(p1, p2);

    Vector rij = p1.getPosition() - p2.getPosition();
    Sim->BCs->applyBC(rij);
    
    return _potential->calculateStepID(rij.nrm() / length_scale);
  }

  double 
  IStepped::getInternalEnergy() const 
  { 
    double Energy = 0.0;
    for (const ICapture::value_type& IDs : *this)
      Energy += getInternalEnergy(Sim->particles[IDs.first.first], Sim->particles[IDs.first.second]);
    return Energy; 
  }

  double 
  IStepped::getInternalEnergy(const Particle& p1, const Particle& p2) const
  {
    ICapture::const_iterator capstat = ICapture::find(ICapture::key_type(p1, p2));
    if (capstat == ICapture::end()) return 0;
    const double energy_scale = _energyScale->getProperty(p1, p2);
    return (*_potential)[capstat->second - 1].second * energy_scale;
  }

  IntEvent
  IStepped::getEvent(const Particle &p1, const Particle &p2) const
  {
  
#ifdef DYNAMO_DEBUG
    if (!Sim->dynamics->isUpToDate(p1))
      M_throw() << "Particle 1 is not up to date";
  
    if (!Sim->dynamics->isUpToDate(p2))
      M_throw() << "Particle 2 is not up to date";

    if (p1 == p2)
      M_throw() << "You shouldn't pass p1==p2 events to the interactions!";
#endif 

    ICapture::const_iterator capstat = ICapture::find(ICapture::key_type(p1, p2));
    const size_t current_step_ID = (capstat == ICapture::end()) ? 0 : capstat->second;
    const std::pair<double, double> step_bounds = _potential->getStepBounds(current_step_ID);
    const double length_scale = _lengthScale->getProperty(p1, p2);

    IntEvent retval(p1, p2, HUGE_VAL, NONE, *this);
    if (step_bounds.first != 0)
      {//Test for the inner step capture
	const double dt = Sim->dynamics->SphereSphereInRoot(p1, p2, step_bounds.first * length_scale);
	if (dt != HUGE_VAL)
	  retval = IntEvent(p1, p2, dt, STEP_IN, *this);
      }

    if (!std::isinf(step_bounds.second))
      {//Test for the outer step capture
	const double dt = Sim->dynamics->SphereSphereOutRoot(p1, p2, step_bounds.second * length_scale);
	if (retval.getdt() > dt)
	  retval = IntEvent(p1, p2, dt, STEP_OUT, *this);
      }
    
    return retval;
  }

  PairEventData
  IStepped::runEvent(Particle& p1, Particle& p2, const IntEvent& iEvent)
  {
    ++Sim->eventCount;

    const double length_scale = _lengthScale->getProperty(p1, p2);
    const double energy_scale = _energyScale->getProperty(p1, p2);

    ICapture::const_iterator capstat = ICapture::find(ICapture::key_type(p1, p2));
    const size_t old_step_ID = (capstat == ICapture::end()) ? 0 : capstat->second;
    const std::pair<double, double> step_bounds = _potential->getStepBounds(old_step_ID);

    size_t new_step_ID;
    size_t edge_ID;
    double diameter;
    switch (iEvent.getType())
      {
      case STEP_OUT:
	{
	  new_step_ID = _potential->outer_step_ID(old_step_ID);
	  edge_ID = _potential->outer_edge_ID(old_step_ID);
	  diameter = step_bounds.second * length_scale;
	  break;
	}
      case STEP_IN:
	{
	  new_step_ID = _potential->inner_step_ID(old_step_ID);
	  edge_ID = _potential->inner_edge_ID(old_step_ID);
	  diameter = step_bounds.first * length_scale;
	  break;
	}
      default:
	M_throw() << "Unknown event type";
      } 

    PairEventData retVal = Sim->dynamics->SphereWellEvent(iEvent, _potential->getEnergyChange(new_step_ID, old_step_ID) * energy_scale, diameter * diameter, new_step_ID);
    EdgeData& data = _edgedata[std::pair<size_t, EEventType>(edge_ID, retVal.getType())];
    ++data.counter;
    data.rdotv_sum += retVal.rvdot;
    //Check if the particles changed their step ID
    if (retVal.getType() != BOUNCE) ICapture::operator[](ICapture::key_type(p1, p2)) = new_step_ID;
    return retVal;
  }

  bool
  IStepped::validateState(const Particle& p1, const Particle& p2, bool textoutput) const
  {
    ICapture::const_iterator capstat = ICapture::find(ICapture::key_type(p1, p2));
    const size_t stored_step_ID = (capstat == ICapture::end()) ? 0 : capstat->second;
    const size_t calculated_step_ID = captureTest(p1, p2);
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
	<< magnet::xml::attr("LengthScale") << _lengthScale->getName()
	<< magnet::xml::attr("EnergyScale") << _energyScale->getName()
	<< *range;
  
    XML << _potential;

    ICapture::outputCaptureMap(XML);  
  }

  void 
  IStepped::outputData(magnet::xml::XmlStream& XML) const
  {
    //This code might not be consistent with units, as the potential
    //is stored in configuration units, not simulation. Also, this
    //code might not be correct if the potential changes direction
    //(id's increase outwards).


    using namespace magnet::xml;
    XML << tag("Interaction")
	<< attr("Name") << intName
	<< attr("Type") << "Stepped"
	<< tag("AccessedSteps")
	<< attr("Direction") << (_potential->direction() ? "Outward" : "Inward")
	<< attr("MaxDiameter") << _potential->max_distance()
      ;

    double kT(0);
    {
      std::shared_ptr<EnsembleNVT> ensemble = std::dynamic_pointer_cast<EnsembleNVT>(Sim->ensemble);
      if (ensemble)
	kT = ensemble->getEnsembleVals()[2];
    }

    for (size_t i(0); i < _potential->cached_steps(); ++i)
      {
	double deltaU = (*_potential)[i].second;
	if (i > 0) deltaU -= (*_potential)[i-1].second;
	
	const double R = (*_potential)[i].first;

	XML << tag("Step") 
	    << attr("ID") << i
	    << attr("R") << R
	    << attr("U") << (*_potential)[i].second
	    << attr("DeltaU") << deltaU
	  ;

	for (const auto& data: _edgedata)
	  if (data.first.first == i)
	    {

	      XML << tag("Event")
		  << attr("Type") << data.first.second
		  << attr("Count") << data.second.counter
		  << attr("RdotV") << data.second.rdotv_sum / (data.second.counter * Sim->units.unitVelocity() * Sim->units.unitLength());
	      
	      if (kT)
		{
		  double gr = 2 * (Sim->getSimVolume() / (4 * R * R * std::sqrt(M_PI * kT) * Sim->N() * Sim->N())) * (data.second.counter / Sim->systemTime);
		  
		  switch (data.first.second)
		    {
		    case EEventType::STEP_OUT:
		      if (deltaU < 0)
			XML << attr("gr") << gr * std::exp(std::abs(deltaU) / kT);
		      else
			XML << attr("gr") << gr;
		      break;
		    case EEventType::STEP_IN:
		      if (deltaU > 0)
			XML << attr("gr") << gr * std::exp(std::abs(deltaU) / kT);
		      else	
			XML << attr("gr") << gr;
		      break;
		    case EEventType::BOUNCE:
		      XML << attr("gr") << - gr * (std::exp(-std::abs(deltaU) / kT) / (std::exp(-std::abs(deltaU) / kT) - 1));
		      break;
		    default:
		      break;
		    }
		}

	      XML << endtag("Event");
	    }
	XML << endtag("Step");
      }
    
    shared_ptr<OPRadialDistribution> raddist = Sim->getOutputPlugin<OPRadialDistribution>();
    if (raddist && kT)
      {
	XML << tag("gr") << chardata();
	const double grBinWidth = raddist->getBinWidth();
	std::vector<std::pair<double, double> > grdata = raddist->getgrdata(0, 0);
	std::vector<std::pair<double, double> > yrdata;
	
	for (size_t i = 1; i < grdata.size(); ++i)
	  {
	    int potential_step = -1;
	    for (size_t stepID(0); stepID < _potential->cached_steps(); ++stepID)
	      {
		const double R = (*_potential)[stepID].first * Sim->units.unitLength();
		const size_t istep = static_cast<size_t>(R / grBinWidth + 0.5);
		if (i == istep) { potential_step = stepID; break; }
	      }
	    
	    if (potential_step == -1)
	      {
		const double R = grdata[i].first / Sim->units.unitLength();
		XML << R << " " << grdata[i].second << "\n";

		if (grdata[i].second == 0)
		  yrdata.push_back(std::pair<double, double>(R, 0));
		else
		  {
		    const size_t step_ID = _potential->calculateStepID(R);
		    double U = 0;
		    if (step_ID > 0) U = (*_potential)[step_ID - 1].second;
		    double yrval = grdata[i].second * std::exp(U * Sim->units.unitEnergy() / kT);
		    if (!std::isnan(yrval))
		      yrdata.push_back(std::pair<double, double>(R, yrval));
		  }
	      }
	    else
	      {
		for (EEventType etype: {EEventType:STEP_OUT, EEventType:BOUNCE, EEventType:STEP_IN})
		  for (const auto& data: _edgedata)
		    if ((data.first.first == potential_step) && (data.first.second == etype))
		      {
			const double R = (*_potential)[potential_step].first;
			double deltaU = (*_potential)[potential_step].second;
			if (potential_step > 0) deltaU -= (*_potential)[potential_step-1].second;
			
			double gr = 2 * (Sim->getSimVolume() / (4 * R * R * std::sqrt(M_PI * kT) * Sim->N() * Sim->N())) * (data.second.counter / Sim->systemTime);
			switch (data.first.second)
			  {
			  case EEventType::STEP_OUT:
			    if (deltaU < 0)
			      gr *= std::exp(std::abs(deltaU) / kT);
			    break;
			  case EEventType::STEP_IN:
			    if (deltaU > 0)
			      gr *= std::exp(std::abs(deltaU) / kT);
			    break;
			  case EEventType::BOUNCE:
			    gr *= -(std::exp(-std::abs(deltaU) / kT) / (std::exp(-std::abs(deltaU) / kT) - 1));
			    break;
			  default: break;
			  }
			XML << R << " " << gr << " *\n";
		      }
	      }
	  }
	XML << endtag("gr")
	    << tag("yr")
	    << chardata();
	for (auto entry : yrdata) XML << entry.first << " " << entry.second << "\n";
	XML << endtag("yr");
      }

    XML << endtag("AccessedSteps")
	<< endtag("Interaction");
    
  }
}
