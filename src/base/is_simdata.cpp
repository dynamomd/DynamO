/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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
    Ensemble(NULL),
    dSysTime(0.0),
    freestreamAcc(0.0),
    lNColl(0),
    lMaxNColl(100000),
    lNPrint(50000),
    lPrintLimiter(0),
    lN(0),
    ptrScheduler(NULL),
    Dynamics(this),
    aspectRatio(1,1,1),
    ranGenerator(static_cast<unsigned>(std::time(0))),
    normal_sampler(ranGenerator, boost::normal_distribution_01<Iflt>()),
    uniform_sampler(ranGenerator),
    lastRunMFT(0.0),
    simID(0),
    status(START),
    binaryXML(true)
  {
#ifdef DYNAMO_CONDOR
    binaryXML = false;
#endif
  }

  SimData::~SimData()
  {
    if (ptrScheduler != NULL) delete ptrScheduler;
  }
  
  void 
  SimData::signalParticleUpdate
  (const CNParticleData& pdat) const
  {
    BOOST_FOREACH(const particleUpdateFunc& func, _particleUpdateNotify)
      func(pdat);
  }

  void 
  SimData::replexerSwap(SimData& other)
  {
    //Get all particles up to date and zero the pecTimes
    Dynamics.Liouvillean().updateAllParticles();
    other.Dynamics.Liouvillean().updateAllParticles();
      
    std::swap(dSysTime, other.dSysTime);
    std::swap(lNColl, other.lNColl);
    std::swap(_particleUpdateNotify, other._particleUpdateNotify);
    
    Dynamics.getSystemEvents().swap(other.Dynamics.getSystemEvents());

    BOOST_FOREACH(smrtPlugPtr<CSystem>& aPtr, Dynamics.getSystemEvents())
      aPtr->changeSystem(this);

    BOOST_FOREACH(smrtPlugPtr<CSystem>& aPtr, other.Dynamics.getSystemEvents())
      aPtr->changeSystem(&other);

    //Rescale the velocities     
    Iflt scale1(sqrt(other.Ensemble->getEnsembleVals()[2]
		     / Ensemble->getEnsembleVals()[2]));
    
    BOOST_FOREACH(CParticle& part, vParticleList)
      part.scaleVelocity(scale1);
    
    other.ptrScheduler->rescaleTimes(scale1);
    
    Iflt scale2(1.0 / scale1);

    BOOST_FOREACH(CParticle& part, other.vParticleList)
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
      std::vector<smrtPlugPtr<COutputPlugin> >::iterator iPtr1 = outputPlugins.begin(), 
	iPtr2 = other.outputPlugins.begin();
      
      while (iPtr1 != outputPlugins.end())
	{
#ifdef DYNAMO_DEBUG
	  if (typeid(*(*iPtr1)) != typeid(*(*iPtr2)))
	    D_throw() << "Output plugin mismatch while replexing! lists not sorted the same perhaps?";
#endif
	  
	  (*iPtr1)->changeSystem(iPtr2->get_ptr());
	  
	  (*iPtr1)->temperatureRescale(scale1);
	  (*iPtr2)->temperatureRescale(scale2);
	  
	  ++iPtr1; 
	  ++iPtr2;
	}
    }

    //This is swapped last as things need it for calcs
    Ensemble->exchange(*other.Ensemble);
  }
}
