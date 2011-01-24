/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
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

#include "sleep.hpp"
#include "globEvent.hpp"
#include "../NparticleEventData.hpp"
#include "../../base/is_simdata.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../../schedulers/scheduler.hpp"

GSleep::GSleep(DYNAMO::SimData* nSim, const std::string& name):
  Global(nSim, "Sleep")
{
  globName = name;
  I_cout() << "Sleep Loaded";
}

GSleep::GSleep(const XMLNode &XML, DYNAMO::SimData* ptrSim):
  Global(ptrSim, "Sleep")
{
  operator<<(XML);

  I_cout() << " Loaded";
}

void 
GSleep::initialise(size_t nID)
{
  ID=nID;

  sleepTime.resize(Sim->N);

  BOOST_FOREACH(const Particle& part, Sim->particleList)
    sleepTime[part.getID()] = HUGE_VAL;

  Sim->registerParticleUpdateFunc
    (magnet::function::MakeDelegate(this, &GSleep::particlesUpdated));
}

void 
GSleep::particlesUpdated(const NEventData& PDat)
{  
  BOOST_FOREACH(const PairEventData& pdat, PDat.L2partChanges)
    {
      const Particle& p1 = pdat.particle1_.getParticle();
      const Particle& p2 = pdat.particle2_.getParticle();
      sleepTime[p1.getID()] = HUGE_VAL;
      sleepTime[p2.getID()] = HUGE_VAL;
    
      //We will assume that there are only two states
      if( p1.testState(Particle::DYNAMIC) !=  p2.testState(Particle::DYNAMIC))
	{
	  //We must assign p1 to the dynamic particle before continuing

	  bool collision = FALSE; //We chech if the event is a collision event
	  //bool convergePos = ; //We chech if the position converges
	  //bool convergeVel = ; //We chech if the velocity converges
	  Vector g(0,0,0);        //We need gravity in order to assure the 
	                          //geometry of the sleeping position  
	  bool Vg = (p1.getVelocity() | g) > 0?TRUE:FALSE; // We need this to be negative, i.e., particle goes down
	  double SleepVel = 1.1 * Sim->dynamics.units().unitVelocity(); //Sleeping velocity, under this you sleep. 
	  //It has to be larger than VElastic
	  double vel = p1.getVelocity().nrm();
	  if(vel < SleepVel)
	    {
	      //std::cerr << "\nRequesting sleep! pID = " 
	      //		<< (p1.testState(Particle::DYNAMIC) ? p1.getID() : p2.getID()) <<"\n";
	      //We sleep just the particle that is awake
	      if (p1.testState(Particle::DYNAMIC))
		sleepTime[p1.getID()] = 0;
	      else
		sleepTime[p2.getID()] = 0;
	    }
	}
    }
}

void 
GSleep::operator<<(const XMLNode& XML)
{
  try {
    globName = XML.getAttribute("Name");	
  }
  catch(...)
    {
      M_throw() << "Error loading GSleep";
    }
}

GlobalEvent 
GSleep::getEvent(const Particle& part) const
{
  return GlobalEvent(part,sleepTime[part.getID()] , VIRTUAL, *this);
}

void 
GSleep::runEvent(const Particle& part) const
{
  GlobalEvent iEvent(getEvent(part));

#ifdef DYNAMO_DEBUG 
  if (boost::math::isnan(iEvent.getdt()))
    M_throw() << "A NAN Interaction collision time has been found"
	      << iEvent.stringData(Sim);
  
  if (iEvent.getdt() == HUGE_VAL)
    M_throw() << "An infinite Interaction (not marked as NONE) collision time has been found\n"
	      << iEvent.stringData(Sim);
#endif

  Sim->dSysTime += iEvent.getdt();
    
  Sim->ptrScheduler->stream(iEvent.getdt());
  
  Sim->dynamics.stream(iEvent.getdt());

  Sim->dynamics.getLiouvillean().updateParticle(part);

  //std::cerr << "\nPerforming sleep! pID = " 
  //	    << part.getID() << "\n";

  //Here is where the particle goes to sleep
  const_cast<Particle&>(part).clearState(Particle::DYNAMIC);
  const_cast<Particle&>(part).getVelocity() = Vector(0,0,0);
  //Sim->particleList[part.getID()].clearState(Particle::DYNAMIC);
  sleepTime[part.getID()] = HUGE_VAL;

  Sim->freestreamAcc += iEvent.getdt();

  Sim->ptrScheduler->fullUpdate(part);
}

void 
GSleep::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Type") << "Sleep"
      << xml::attr("Name") << globName;
}
