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
namespace {
  class densityCalculator
  {
  public:
    densityCalculator() { density = 0; }
    
    void foundparticle(Particle p1, size_t p2ID)
    { 
      density += 1;	
    } 
    
    double density;
  };

}

double GSleep::getDensity(const Particle& part)
{
  densityCalculator dcalc;

  
  // nblist.getParticleNeighbourhood(particle, &densityCalculator::foundparticle, &dcalc);
  return dcalc.density;
}

GSleep::GSleep(DYNAMO::SimData* nSim,  CRange* range, const std::string& name,const double sl):
  Global(range, nSim, "Sleep"),
  sleepVelocity(sl)
{
  globName = name;
  I_cout() << "Sleep Loaded";
}

GSleep::GSleep(const XMLNode &XML, DYNAMO::SimData* ptrSim):
  Global(NULL, ptrSim, "Sleep")
{
  operator<<(XML);

  I_cout() << " Loaded";
}

void 
GSleep::initialise(size_t nID)
{

  //Neighbour
  // try {
  //   NBListID = Sim->dynamics.getGlobal("SchedulerNBList")->getID();
  // }
  // catch(std::exception& cxp)
  //   {
  //     M_throw() << "Failed while finding the neighbour list global.\n"
  // 		<< "You must have a neighbour list enabled for this\n"
  // 		<< "scheduler called SchedulerNBList.\n"
  // 		<< cxp.what();
  //   }
  
  // if (dynamic_cast<const CGNeighbourList*>
  //     (Sim->dynamics.getGlobals()[NBListID].get_ptr())
  //     == NULL)
  //   M_throw() << "The Global named SchedulerNBList is not a neighbour list!";

  // static_cast<CGNeighbourList&>
  //   (*Sim->dynamics.getGlobals()[NBListID].get_ptr())
  //   .markAsUsedInScheduler();
  //end neighbour creation

  ID=nID;

  lastPosition.resize(Sim->N);
  lastVelocity.resize(Sim->N);
  Vector aux(0,0,0);
  BOOST_FOREACH(const Particle& part, Sim->particleList)
    {
      lastVelocity[part.getID()] = aux;
      lastPosition[part.getID()] = aux;
    }

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
    
      //We will assume that there are only two states
      if (range->isInRange(p1) || range->isInRange(p2))
	if (p1.testState(Particle::DYNAMIC) !=  p2.testState(Particle::DYNAMIC))
	{
	  //This is the particle which is dynamic
	  const Particle& dp = p1.testState(Particle::DYNAMIC) ? p1 : p2;
	  const Particle& sp = p1.testState(Particle::DYNAMIC) ? p2 : p1;
	  	  
	  double vel = dp.getVelocity().nrm();

	  //If the static particle is in range, resleep it and mark its momentum to be transferred to the other particle 
	  if (range->isInRange(sp))
	    {
	      stateChange[sp.getID()] = zeroedVector();
	      stateChange[dp.getID()] -= sp.getVelocity() * Sim->dynamics.getSpecies(sp).getMass();
	    }
	  
	  if (range->isInRange(dp))
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
	      if((vel < sleepVelocity) &&  Vg && convergeVel && convergePos)
		stateChange[dp.getID()] = zeroedVector();
	    }

	  double normalVel = dp.getVelocity()|(dp.getPosition()-sp.getPosition());
	    
	  lastVelocity[p1.getID()] = p1.getVelocity();
	  lastVelocity[p2.getID()] = p2.getVelocity();
	  lastPosition[p1.getID()] = p1.getPosition();
	  lastPosition[p2.getID()] = p2.getPosition();
	}
    }
}

void 
GSleep::operator<<(const XMLNode& XML)
{
  range.set_ptr(CRange::loadClass(XML, Sim));

  try {
    globName = XML.getAttribute("Name");
    sleepVelocity = Sim->dynamics.units().unitVelocity() * 
      boost::lexical_cast<double>(XML.getAttribute("SleepV"));
  }
  catch(...)
    {
      M_throw() << "Error loading GSleep";
    }
}

GlobalEvent
GSleep::getEvent(const Particle& part) const
{
  if (stateChange.find(part.getID()) != stateChange.end())//Check if we want a state change
    return GlobalEvent(part, 0, CHECK, *this);

  return GlobalEvent(part, HUGE_VAL, NONE, *this);
}

void 
GSleep::runEvent(const Particle& part, const double dt) const
{
  GlobalEvent iEvent(getEvent(part));
  iEvent.setdt(dt); //We only trust the schedulers time, as we don't
		    //track the motion of the system in Globals
  
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

  //Here is where the particle goes to sleep or wakes
  ++Sim->eventCount;
  ParticleEventData EDat(part, Sim->dynamics.getSpecies(part), iEvent.getType());

  if (part.testState(Particle::DYNAMIC))
    if (stateChange[part.getID()].nrm() == 0)
      {
	const_cast<Particle&>(part).clearState(Particle::DYNAMIC);
	const_cast<Particle&>(part).getVelocity() = Vector(0,0,0);
	iEvent.setType(SLEEP);
      }
    else
      {
	const_cast<Particle&>(part).getVelocity() += stateChange[part.getID()] / EDat.getSpecies().getMass();
	iEvent.setType(CHECK);
      }
  else
    {
      const_cast<Particle&>(part).getVelocity() = Vector(0,0,0);
      iEvent.setType(RESLEEP);
    }
  stateChange.erase(part.getID());

  EDat.setDeltaKE(0.5 * EDat.getSpecies().getMass()
		  * (part.getVelocity().nrm2() 
		     - EDat.getOldVel().nrm2()));

  Sim->signalParticleUpdate(EDat);

  //Now we're past the event, update the scheduler and plugins
  Sim->ptrScheduler->fullUpdate(part);

  BOOST_FOREACH(magnet::ClonePtr<OutputPlugin> & Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(iEvent, EDat);

}

void 
GSleep::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Type") << "Sleep"
      << xml::attr("Name") << globName
      << xml::attr("SleepV") << sleepVelocity/Sim->dynamics.units().unitVelocity()
      << range;
}
