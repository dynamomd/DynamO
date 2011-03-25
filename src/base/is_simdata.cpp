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

#include "is_simdata.hpp"
#include "../schedulers/scheduler.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"
#include "../../schedulers/scheduler.hpp"
#include "../dynamics/systems/system.hpp"

namespace DYNAMO
{
  SimData::SimData():
    ensemble(NULL),
    dSysTime(0.0),
    freestreamAcc(0.0),
    eventCount(0),
    endEventCount(100000),
    eventPrintInterval(50000),
    nextPrintEvent(0),
    N(0),
    ptrScheduler(NULL),
    dynamics(this),
    aspectRatio(1,1,1),
    ranGenerator(static_cast<unsigned>(std::time(0))),
    normal_sampler(ranGenerator, boost::normal_distribution_01<double>()),
    uniform_sampler(ranGenerator),
    lastRunMFT(0.0),
    simID(0),
    replexExchangeNumber(0),
    status(START),
    binaryXML(true)
  {
  }

  SimData::~SimData()
  {
    if (ptrScheduler != NULL) delete ptrScheduler;
  }
  
  void 
  SimData::signalParticleUpdate
  (const NEventData& pdat) const
  {
    BOOST_FOREACH(const particleUpdateFunc& func, _particleUpdateNotify)
      func(pdat);
  }

  void 
  SimData::replexerSwap(SimData& other)
  {
    //Get all particles up to date and zero the pecTimes
    dynamics.getLiouvillean().updateAllParticles();
    other.dynamics.getLiouvillean().updateAllParticles();
      
    std::swap(dSysTime, other.dSysTime);
    std::swap(eventCount, other.eventCount);
    std::swap(_particleUpdateNotify, other._particleUpdateNotify);
    
    dynamics.getSystemEvents().swap(other.dynamics.getSystemEvents());

    BOOST_FOREACH(magnet::ClonePtr<System>& aPtr, dynamics.getSystemEvents())
      aPtr->changeSystem(this);

    BOOST_FOREACH(magnet::ClonePtr<System>& aPtr, other.dynamics.getSystemEvents())
      aPtr->changeSystem(&other);

    //Rescale the velocities 
    double scale1(sqrt(other.ensemble->getEnsembleVals()[2]
		     / ensemble->getEnsembleVals()[2]));
    
    BOOST_FOREACH(Particle& part, particleList)
      part.scaleVelocity(scale1);
    
    other.ptrScheduler->rescaleTimes(scale1);
    
    double scale2(1.0 / scale1);

    BOOST_FOREACH(Particle& part, other.particleList)
      part.scaleVelocity(scale2);
    
    ptrScheduler->rescaleTimes(scale2);
    
    ptrScheduler->rebuildSystemEvents();
    other.ptrScheduler->rebuildSystemEvents();    

    //Globals?
#ifdef DYNAMO_DEBUG
    if (outputPlugins.size() != other.outputPlugins.size())
      std::cerr << "Error, could not swap output plugin lists as they are not equal in size";
#endif

    outputPlugins.swap(other.outputPlugins);      
    
    {
      std::vector<magnet::ClonePtr<OutputPlugin> >::iterator iPtr1 = outputPlugins.begin(), 
	iPtr2 = other.outputPlugins.begin();
      
      while (iPtr1 != outputPlugins.end())
	{
#ifdef DYNAMO_DEBUG
	  if (typeid(*(*iPtr1)) != typeid(*(*iPtr2)))
	    M_throw() << "Output plugin mismatch while replexing! lists not sorted the same perhaps?";
#endif
	  
	  (*iPtr1)->changeSystem(iPtr2->get_ptr());
	  
	  (*iPtr1)->temperatureRescale(scale1 * scale1);
	  (*iPtr2)->temperatureRescale(scale2 * scale2);
	  
	  ++iPtr1; 
	  ++iPtr2;
	}
    }

    //This is swapped last as things need it for calcs
    ensemble->exchange(*other.ensemble);
  }
}
