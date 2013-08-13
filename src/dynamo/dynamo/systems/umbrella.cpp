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

#include <dynamo/systems/umbrella.hpp>

#include <dynamo/units/units.hpp>
#include <dynamo/BC/BC.hpp>
#include <dynamo/particle.hpp>
#include <dynamo/species/species.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/ranges/include.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/lexical_cast.hpp>

#ifdef DYNAMO_DEBUG 
#include <boost/math/special_functions/fpclassify.hpp>
#endif

namespace dynamo {

  SysUmbrella::SysUmbrella(const magnet::xml::Node& XML, dynamo::Simulation* tmp): 
    System(tmp),
    _stepID(std::numeric_limits<size_t>::max())
  {
    dt = HUGE_VAL;
    operator<<(XML);
    type = UMBRELLA;
  }

  void 
  SysUmbrella::runEvent() const
  {
    double locdt = dt;
  
#ifdef DYNAMO_DEBUG 
    if (boost::math::isnan(locdt))
      M_throw() << "A NAN system event time has been found";
#endif
  
    Sim->systemTime += locdt;
  
    Sim->ptrScheduler->stream(locdt);
  
    //dynamics must be updated first
    Sim->stream(locdt);

    ++Sim->eventCount;

    for (const size_t& id : *range1)
      Sim->dynamics->updateParticle(Sim->particles[id]);
  
    for (const size_t& id : *range2)
      Sim->dynamics->updateParticle(Sim->particles[id]);

    size_t new_step_ID;
    if (type == STEP_OUT)
      new_step_ID = _stepID - 1 + 2 * _potential->direction();
    else if (type == STEP_IN)
      new_step_ID = _stepID + 1 - 2 * _potential->direction();
    else
      M_throw() << "Unknown event type";
  
    EEventType etype(NONE);
    NEventData SDat(Sim->dynamics->multibdyWellEvent(*range1, *range2, 0.0, _potential->getEnergyChange(new_step_ID, _stepID) * _energyScale, etype));

    if (etype != BOUNCE) _stepID = new_step_ID;

    Sim->_sigParticleUpdate(SDat);
  
    //Only 1ParticleEvents occur
    for (const ParticleEventData& PDat : SDat.L1partChanges)
      Sim->ptrScheduler->fullUpdate(Sim->particles[PDat.getParticleID()]);
  
    for (shared_ptr<OutputPlugin>& Ptr : Sim->outputPlugins)
      Ptr->eventUpdate(*this, SDat, locdt); 
  }

  void
  SysUmbrella::initialise(size_t nID)
  {
    ID = nID;

    if (_stepID == std::numeric_limits<size_t>::max())
      {
	for(const size_t& id : *range1)
	  Sim->dynamics->updateParticle(Sim->particles[id]);
	
	for(const size_t& id : *range2)
	  Sim->dynamics->updateParticle(Sim->particles[id]);
	
	std::pair<Vector, Vector> r1data = Sim->dynamics->getCOMPosVel(*range1);
	std::pair<Vector, Vector> r2data = Sim->dynamics->getCOMPosVel(*range2);
	Vector r12 = r1data.first - r2data.first;
	Sim->BCs->applyBC(r12);
	_stepID = _potential->calculateStepID(r12.nrm() / _lengthScale);
      }
  
    recalculateTime();

    Sim->_sigParticleUpdate.connect<SysUmbrella, &SysUmbrella::particlesUpdated>(this);
  }

  void 
  SysUmbrella::recalculateTime()
  {
    for (const size_t& id : *range1)
      Sim->dynamics->updateParticle(Sim->particles[id]);
  
    for (const size_t& id : *range2)
      Sim->dynamics->updateParticle(Sim->particles[id]);
  
    dt = HUGE_VAL;
    type = NONE;

    const std::pair<double, double> step_bounds = _potential->getStepBounds(_stepID);
    
    if (step_bounds.first != 0)
      {
	const double new_dt = Sim->dynamics->SphereSphereInRoot(*range1, *range2, step_bounds.first * _lengthScale);
	if (new_dt < dt)
	  {
	    dt = new_dt;
	    type = STEP_IN;
	  }
      }
    
    if (!std::isinf(step_bounds.second))
      {
	const double new_dt = Sim->dynamics->SphereSphereOutRoot(*range1, *range2, step_bounds.second * _lengthScale);

	if (new_dt < dt)
	  {
	    dt = new_dt;
	    type = STEP_OUT;	    
	  }
      }
  }

  void 
  SysUmbrella::particlesUpdated(const NEventData& PDat)
  {
    for (const ParticleEventData& pdat : PDat.L1partChanges)
      if (range1->isInRange(Sim->particles[pdat.getParticleID()])
	  || range2->isInRange(Sim->particles[pdat.getParticleID()]))
	{
	  recalculateTime();
	  Sim->ptrScheduler->rebuildSystemEvents();
	  return;
	}

    for (const PairEventData& pdat : PDat.L2partChanges)
      if (range1->isInRange(Sim->particles[pdat.particle1_.getParticleID()])
	  || range2->isInRange(Sim->particles[pdat.particle1_.getParticleID()])
	  || range1->isInRange(Sim->particles[pdat.particle2_.getParticleID()])
	  || range2->isInRange(Sim->particles[pdat.particle2_.getParticleID()]))
	{
	  recalculateTime();
	  Sim->ptrScheduler->rebuildSystemEvents();
	  return;
	}
  }

  void
  SysUmbrella::operator<<(const magnet::xml::Node& XML)
  {
    sysName = XML.getAttribute("Name");
    magnet::xml::Node rangeNode = XML.getNode("IDRange");
    range1 = shared_ptr<IDRange>(IDRange::getClass(rangeNode, Sim));
    ++rangeNode;
    range2 = shared_ptr<IDRange>(IDRange::getClass(rangeNode, Sim));
    _potential = Potential::getClass(XML.getNode("Potential"));

    _lengthScale = XML.getAttribute("LengthScale").as<double>() * Sim->units.unitLength();
    _energyScale = XML.getAttribute("EnergyScale").as<double>() * Sim->units.unitEnergy();

    if (XML.hasAttribute("CurrentStep"))
      _stepID = XML.getAttribute("CurrentStep").as<size_t>();
  }

  void 
  SysUmbrella::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::tag("System")
	<< magnet::xml::attr("Type") << "Umbrella"
	<< magnet::xml::attr("Name") << sysName
	<< magnet::xml::attr("LengthScale") << _lengthScale / Sim->units.unitLength()
	<< magnet::xml::attr("EnergyScale") << _energyScale / Sim->units.unitEnergy()
	<< magnet::xml::attr("CurrentStep") << _stepID
	<< _potential
	<< range1
	<< range2
	<< magnet::xml::endtag("System");
  }
}
