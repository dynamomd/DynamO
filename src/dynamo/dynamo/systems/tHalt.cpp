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

#include <dynamo/systems/tHalt.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>

namespace dynamo {
  SystHalt::SystHalt(dynamo::Simulation* nSim, double ndt, std::string nName):
    System(nSim)
  {
    dt = ndt * Sim->units.unitTime();

    sysName = nName;

    dout << "System halt set for " 
	 << ndt << std::endl;

    type = VIRTUAL;
  }

  void
  SystHalt::runEvent()
  {
    Event event = getEvent();
  
#ifdef DYNAMO_DEBUG 
    if (std::isnan(event._dt))
      M_throw() << "A NAN system event time has been found";
#endif
    
    Sim->systemTime += event._dt;
    
    Sim->ptrScheduler->stream(event._dt);
  
    //dynamics must be updated first
    Sim->stream(event._dt);

    NEventData SDat;
    Sim->_sigParticleUpdate(SDat);
    
    for (shared_ptr<OutputPlugin>& Ptr : Sim->outputPlugins)
      Ptr->eventUpdate(event, SDat);
  
    Sim->nextPrintEvent = Sim->endEventCount = Sim->eventCount;
  }

  void 
  SystHalt::initialise(size_t nID)
  {
    ID=nID;
  }

  void 
  SystHalt::setdt(double ndt)
  { dt = ndt * Sim->units.unitTime(); }

  void 
  SystHalt::increasedt(double ndt)
  { 
    dt += ndt * Sim->units.unitTime(); 
  }
}
