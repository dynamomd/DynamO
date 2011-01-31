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

#include "sleep.hpp"
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random/uniform_int.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../dynamics.hpp"
#include "../units/units.hpp"
#include "../BC/BC.hpp"
#include "../../simulation/particle.hpp"
#include "../species/species.hpp"
#include "../NparticleEventData.hpp"
#include "../ranges/include.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../../schedulers/scheduler.hpp"

SSleep::SSleep(const XMLNode& XML, DYNAMO::SimData* tmp): 
  System(tmp),
  _range(NULL)
{
  dt = HUGE_VAL;
  operator<<(XML);
  type = SLEEP;
}

SSleep::SSleep(DYNAMO::SimData* nSim, std::string nName, CRange* r1, double sleepV):
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

  lastPosition.resize(Sim->N);
  lastVelocity.resize(Sim->N);
  Vector aux(0,0,0);
  BOOST_FOREACH(const Particle& part, Sim->particleList)
    {
      lastVelocity[part.getID()] = aux;
      lastPosition[part.getID()] = aux;
    }

  recalculateTime();
}

void
SSleep::operator<<(const XMLNode& XML)
{
  if (strcmp(XML.getAttribute("Type"),"Sleep"))
    M_throw() << "Attempting to load Sleep from a " 
	      << XML.getAttribute("Type") <<  " entry"; 
  
  try {
    sysName = XML.getAttribute("Name");
    
    _sleepVelocity = Sim->dynamics.units().unitVelocity() * 
      boost::lexical_cast<double>(XML.getAttribute("SleepV"));

    _range.set_ptr(CRange::loadClass(XML, Sim));
  }
  catch (boost::bad_lexical_cast &)
    { M_throw() << "Failed a lexical cast in SSleep"; }
}

void 
SSleep::outputXML(xml::XmlStream& XML) const
{
  XML << xml::tag("System")
      << xml::attr("Type") << "Sleep"
      << xml::attr("Name") << sysName
      << xml::attr("SleepV") << _sleepVelocity / Sim->dynamics.units().unitVelocity()
      << _range
      << xml::endtag("System");
}


void 
SSleep::recalculateTime()
{
  dt = (stateChange.empty()) ? HUGE_VAL : 0;
  type = (stateChange.empty()) ? NONE : SLEEP;
}

void 
SSleep::particlesUpdated(const NEventData& PDat)
{
  BOOST_FOREACH(const PairEventData& pdat, PDat.L2partChanges)
    {
      const Particle& p1 = pdat.particle1_.getParticle();
      const Particle& p2 = pdat.particle2_.getParticle();
    
      //We will assume that there are only two states
      if (_range->isInRange(p1) || _range->isInRange(p2))
	if (p1.testState(Particle::DYNAMIC) !=  p2.testState(Particle::DYNAMIC))
	{
	  //This is the particle which is dynamic
	  const Particle& dp = p1.testState(Particle::DYNAMIC) ? p1 : p2;
	  const Particle& sp = p1.testState(Particle::DYNAMIC) ? p2 : p1;
	  	  
	  double vel = dp.getVelocity().nrm();

	  //If the static particle is in range, resleep it and mark its momentum to be transferred to the other particle 
	  if (_range->isInRange(sp))
	    {
	      stateChange[sp.getID()] = zeroedVector();
	      stateChange[dp.getID()] -= sp.getVelocity() * Sim->dynamics.getSpecies(sp).getMass();
	    }
	  
	  if (_range->isInRange(dp))
	    {
	      bool collision = FALSE; //We chech if the event is a collision event
	      bool convergePos = FALSE; //We chech if the position converges
	      bool convergeVel = FALSE; //We chech if the velocity converges
	      
	      Vector g(0,0,-1);        //We need gravity in order to assure the 
	      //geometry of the sleeping position 
	      
	      //It has to be larger than ElasticV, it needs to be added from command line
	      double converge = 0.01;
	      double wakeUpVel = 0.1;
	      // Here we check the last velocity MARCUS!!!
	      double aux = (dp.getVelocity() - lastVelocity[dp.getID()])|g;
	      if (aux < converge && aux > 0) // Small and converging (>0)
		convergeVel = TRUE;
	      
	      //Position
	      if(((dp.getPosition()-lastPosition[dp.getID()])|g) < converge)
		convergePos = TRUE;
	      
	      // We need this to be negative, i.e., particle goes down
	      bool Vg = (dp.getVelocity() | g) > 0;

	      //If the dynamic particle is going to fall asleep, mark its impulse as 0
	      if((vel < _sleepVelocity) &&  Vg && convergeVel && convergePos)
		stateChange[dp.getID()] = zeroedVector();
	    }

	  double normalVel = dp.getVelocity()|(dp.getPosition()-sp.getPosition());
	    
	  lastVelocity[p1.getID()] = p1.getVelocity();
	  lastVelocity[p2.getID()] = p2.getVelocity();
	  lastPosition[p1.getID()] = p1.getPosition();
	  lastPosition[p2.getID()] = p2.getPosition();
	}
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
  double locdt = dt;

  dt = HUGE_VAL;
    
#ifdef DYNAMO_DEBUG 
  if (boost::math::isnan(locdt))
    M_throw() << "A NAN system event time has been found";
#endif
  
  Sim->dSysTime += locdt;
  Sim->ptrScheduler->stream(locdt);
  
  //dynamics must be updated first
  Sim->dynamics.stream(locdt);

  ++Sim->eventCount;

  NEventData SDat;

  typedef std::map<size_t, zeroedVector>::value_type locPair;
  BOOST_FOREACH(const locPair& p, stateChange)
    {
      const Particle& part = Sim->particleList[p.first];
      Sim->dynamics.getLiouvillean().updateParticle(part);
      
      ParticleEventData EDat(part, Sim->dynamics.getSpecies(part), SLEEP);

      if (part.testState(Particle::DYNAMIC))
	if (stateChange[part.getID()].nrm() == 0)
	  {
	    const_cast<Particle&>(part).clearState(Particle::DYNAMIC);
	    const_cast<Particle&>(part).getVelocity() = Vector(0,0,0);
	  }
	else
	  const_cast<Particle&>(part).getVelocity() += stateChange[part.getID()] / EDat.getSpecies().getMass();
      else
	const_cast<Particle&>(part).getVelocity() = Vector(0,0,0);

      EDat.setDeltaKE(0.5 * EDat.getSpecies().getMass()
		      * (part.getVelocity().nrm2() 
			 - EDat.getOldVel().nrm2()));

      SDat.L1partChanges.push_back(EDat);
    }

  Sim->signalParticleUpdate(SDat);

  BOOST_FOREACH(const ParticleEventData& PDat, SDat.L1partChanges)
    Sim->ptrScheduler->fullUpdate(PDat.getParticle());
  
  stateChange.clear();
  
  locdt += Sim->freestreamAcc;

  Sim->freestreamAcc = 0;
  
  BOOST_FOREACH(magnet::ClonePtr<OutputPlugin>& Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(*this, SDat, locdt); 
}

