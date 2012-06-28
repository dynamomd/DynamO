/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>
    Copyright (C) 2011  Sebastian Gonzalez <tsuresuregusa@gmail.com>

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

#include <dynamo/systems/sleep.hpp>

#include <dynamo/units/units.hpp>
#include <dynamo/BC/BC.hpp>
#include <dynamo/particle.hpp>
#include <dynamo/species/species.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/ranges/include.hpp>
#include <dynamo/dynamics/gravity.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random/uniform_int.hpp>

#ifdef DYNAMO_DEBUG 
#include <boost/math/special_functions/fpclassify.hpp>
#endif

namespace dynamo {
  SSleep::SSleep(const magnet::xml::Node& XML, dynamo::Simulation* tmp): 
    System(tmp)
  {
    dt = HUGE_VAL;
    operator<<(XML);
    type = SLEEP;
  }

  SSleep::SSleep(dynamo::Simulation* nSim, std::string nName, Range* r1, double sleepV):
    System(nSim),
    _range(r1),
    _sleepVelocity(sleepV)
  {
    sysName = nName;
    type = SLEEP;
  }

  void
  SSleep::initialise(size_t nID)
  {
    ID = nID;

    Sim->registerParticleUpdateFunc
      (magnet::function::MakeDelegate(this, &SSleep::particlesUpdated));

    recalculateTime();

    _lastData.resize(Sim->N);
    BOOST_FOREACH(const Particle& part, Sim->particleList)
      {
	_lastData[part.getID()].first = part.getPosition();
	_lastData[part.getID()].second = - HUGE_VAL;
      }
  }

  void
  SSleep::operator<<(const magnet::xml::Node& XML)
  {
    if (strcmp(XML.getAttribute("Type"), "Sleep"))
      M_throw() << "Attempting to load Sleep from a " 
		<< XML.getAttribute("Type") <<  " entry"; 
  
    try {
      sysName = XML.getAttribute("Name");
    
      _sleepVelocity = XML.getAttribute("SleepV").as<double>() * Sim->units.unitVelocity();
      _sleepDistance = Sim->units.unitLength() * 0.01;
      _sleepTime = Sim->units.unitTime() * 0.0001;
      _range = shared_ptr<Range>(Range::getClass(XML, Sim));
    }
    catch (boost::bad_lexical_cast &)
      { M_throw() << "Failed a lexical cast in SSleep"; }
  }

  void 
  SSleep::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::tag("System")
	<< magnet::xml::attr("Type") << "Sleep"
	<< magnet::xml::attr("Name") << sysName
	<< magnet::xml::attr("SleepV") << _sleepVelocity / Sim->units.unitVelocity()
	<< _range
	<< magnet::xml::endtag("System");
  }


  void 
  SSleep::recalculateTime()
  {
    dt = (stateChange.empty()) ? HUGE_VAL : -std::numeric_limits<float>::max();
    type = (stateChange.empty()) ? NONE : SLEEP;
  }

  bool 
  SSleep::sleepCondition(const Particle& part, const Vector& g, const Vector& vel)
  {
    Vector diff(part.getPosition() - _lastData[part.getID()].first);
    Sim->BCs->applyBC(diff);
  
    double gnrm = g.nrm();

    return (gnrm != 0)
      && (diff.nrm() < _sleepDistance)
      && ((Sim->dSysTime - _lastData[part.getID()].second) < _sleepTime)
      && (((part.getVelocity() + vel) | (g / gnrm)) < _sleepVelocity);
  }

  void 
  SSleep::particlesUpdated(const NEventData& PDat)
  {
    BOOST_FOREACH(const PairEventData& pdat, PDat.L2partChanges)
      {
	const Particle& p1 = Sim->particleList[pdat.particle1_.getParticleID()];
	const Particle& p2 = Sim->particleList[pdat.particle2_.getParticleID()];
      
	//FC = Fixed collider, DP = Dynamic particle, SP = Static
	//particle, ODP = other dynamic particle, OSP = other static
	//particle

	//Check if either particle is static and should not move

	//[O?P-O?P]
	if (!_range->isInRange(p1) && !_range->isInRange(p2)) continue;

	//DP-[DP/ODP]
	if (p1.testState(Particle::DYNAMIC) && p2.testState(Particle::DYNAMIC)) continue; 
      
	//SP-[FC/SP/OSP]
#ifdef DYNAMO_DEBUG
	if (!p1.testState(Particle::DYNAMIC) && !p2.testState(Particle::DYNAMIC))
	  M_throw() << "Static particles colliding!";
#endif

	//We're guaranteed by the previous tests that
	//p1.testState(Particle::DYNAMIC) != p2.testState(Particle::DYNAMIC)
	//and that at least one particle is in the range

	//Sort the particles
	const Particle& dp = p1.testState(Particle::DYNAMIC) ? p1 : p2;
	const Particle& sp = p1.testState(Particle::DYNAMIC) ? p2 : p1;

	//Other variables
	Vector g(static_cast<const DynGravity&>(*Sim->dynamics).getGravityVector());

	if (!_range->isInRange(sp))  //DP-FC
	  {
	    //If the dynamic particle is going to fall asleep, mark its impulse as 0
	    if (sleepCondition(dp, g))
	      stateChange[dp.getID()] = Vector(0,0,0);
	    continue;
	  }

	if (!_range->isInRange(dp)) continue;

	//Final case
	//DP-SP
	//sp is in the range (a wakeable particle)

	//If the static particle sleeps
	if ((sleepCondition(sp, g)))
	  {
	    double massRatio = Sim->species[sp]->getMass(sp.getID()) 
	      / Sim->species[dp]->getMass(dp.getID());

	    stateChange[sp.getID()] = Vector(0,0,0);
	    stateChange[dp.getID()] = -sp.getVelocity() * massRatio;
	  
	    //Check if the sleep conditions match
	    if ((sleepCondition(dp, g, -sp.getVelocity() * massRatio)))
	      {
		stateChange[dp.getID()] = Vector(0,0,0);
		continue;
	      }

	    //We're not going to sleep the particle due to the standard
	    //rule, but here we try to catch anything that might cause a
	    //problem.

	    //Sometimes our relative velocity effectively goes to zero
	    //(in comparison to the other components). This means the
	    //particle will just keep having an event, we sleep it
	    //instead.
	    if ((pdat.dP.nrm() / Sim->species[dp]->getMass(dp.getID())) 
		< _sleepVelocity)
	      {
		stateChange[dp.getID()] = Vector(0,0,0);
		continue;
	      }
	    
	    continue;
	  }

	//Finally, just wake up the static particle
	stateChange[sp.getID()] = Vector(1,1,1);	
      }

    BOOST_FOREACH(const PairEventData& pdat, PDat.L2partChanges)
      {
	const size_t& p1 = pdat.particle1_.getParticleID();
	_lastData[p1].first = Sim->particleList[p1].getPosition();
	_lastData[p1].second = Sim->dSysTime;

	const size_t& p2 = pdat.particle2_.getParticleID();
	_lastData[p2].first = Sim->particleList[p2].getPosition();
	_lastData[p2].second = Sim->dSysTime;
      }

    if (!stateChange.empty())
      {
	recalculateTime();
	Sim->ptrScheduler->rebuildSystemEvents();
      }
  }

  void 
  SSleep::runEvent() const
  {
    double locdt = 0;

    dt = HUGE_VAL;
    
#ifdef DYNAMO_DEBUG 
    if (boost::math::isnan(locdt))
      M_throw() << "A NAN system event time has been found";
#endif
  
    Sim->dSysTime += locdt;
    Sim->ptrScheduler->stream(locdt);
  
    //dynamics must be updated first
    Sim->stream(locdt);

    //++Sim->eventCount;

    NEventData SDat;

    typedef std::map<size_t, Vector>::value_type locPair;
    BOOST_FOREACH(const locPair& p, stateChange)
      {
	Particle& part = Sim->particleList[p.first];
	Sim->dynamics->updateParticle(part);
      
#ifdef DYNAMO_DEBUG 
	if (stateChange.find(part.getID()) == stateChange.end())
	  M_throw() << "Running an event for a particle with no state change!";
#endif

	EEventType type = WAKEUP;
	if ((stateChange[part.getID()][0] == 0) 
	    && (stateChange[part.getID()][1] == 0) 
	    && (stateChange[part.getID()][2] == 0))
	  {
	    if (part.testState(Particle::DYNAMIC)) 
	      type = SLEEP;
	    else
	      type = RESLEEP;
	  }
	else
	  {
	    if (part.testState(Particle::DYNAMIC)) 
	      type = CORRECT;
	    else
	      type = WAKEUP;
	  }

	ParticleEventData EDat(part, *Sim->species[part], type);

	switch (type)
	  {
	  case SLEEP:
	    part.clearState(Particle::DYNAMIC);
	  case RESLEEP:
	    part.getVelocity() = Vector(0,0,0);
	    break;
	  case CORRECT:
	    part.getVelocity() += stateChange[part.getID()];
	  case WAKEUP:
	    part.setState(Particle::DYNAMIC);
	    break;
	  default:
	    M_throw() << "Bad event type!";
	  }
	  
	EDat.setDeltaKE(0.5 * Sim->species[EDat.getSpeciesID()]->getMass(part.getID())
			* (part.getVelocity().nrm2() 
			   - EDat.getOldVel().nrm2()));

	SDat.L1partChanges.push_back(EDat);
      }

    //Must clear the state before calling the signal, otherwise this
    //will erroneously schedule itself again
    stateChange.clear(); 
    Sim->signalParticleUpdate(SDat);

    BOOST_FOREACH(const ParticleEventData& PDat, SDat.L1partChanges)
      Sim->ptrScheduler->fullUpdate(Sim->particleList[PDat.getParticleID()]);
  
    locdt += Sim->freestreamAcc;

    Sim->freestreamAcc = 0;
  
    BOOST_FOREACH(shared_ptr<OutputPlugin>& Ptr, Sim->outputPlugins)
      Ptr->eventUpdate(*this, SDat, locdt); 
  }
}
