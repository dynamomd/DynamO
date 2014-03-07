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

#include <dynamo/systems/sysTicker.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/outputplugins/tickerproperty/ticker.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/schedulers/scheduler.hpp>

namespace dynamo {
  SysTicker::SysTicker(dynamo::Simulation* nSim, double nPeriod, std::string nName):
    System(nSim)
  {
    if (nPeriod <= 0.0)
      nPeriod = Sim->units.unitTime();

    dt = nPeriod;
    period = nPeriod;

    sysName = nName;

    dout << "System ticker set for a peroid of " 
	 << nPeriod / Sim->units.unitTime() << std::endl;
  }

  void
  SysTicker::runEvent()
  {
    double locdt = dt;
  
#ifdef DYNAMO_DEBUG 
    if (std::isnan(dt))
      M_throw() << "A NAN system event time has been found";
#endif
    

    Sim->systemTime += locdt;

    Sim->ptrScheduler->stream(locdt);
  
    //dynamics must be updated first
    Sim->stream(locdt);
  
    dt += period;
  
    //This is done here as most ticker properties require it
    Sim->dynamics->updateAllParticles();

    for (shared_ptr<OutputPlugin>& Ptr : Sim->outputPlugins)
      {
	shared_ptr<OPTicker> ptr = std::dynamic_pointer_cast<OPTicker>(Ptr);
	if (ptr) ptr->ticker();
      }

    for (shared_ptr<OutputPlugin>& Ptr : Sim->outputPlugins)
      Ptr->eventUpdate(*this, NEventData(), locdt);
  }

  void 
  SysTicker::initialise(size_t nID)
  { ID = nID; }

  void 
  SysTicker::setdt(double ndt)
  { 
    dt = ndt * Sim->units.unitTime(); 
  }

  void 
  SysTicker::increasedt(double ndt)
  { 
    dt += ndt * Sim->units.unitTime(); 
  }

  void 
  SysTicker::setTickerPeriod(const double& nP)
  { 
    dout << "Setting system ticker period to " 
	 << nP / Sim->units.unitTime() << std::endl;

    period = nP; 

    dt = nP;

    if ((Sim->status >= INITIALISED) && Sim->endEventCount)
      Sim->ptrScheduler->rebuildSystemEvents();
  }
}
