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

#include <dynamo/systems/schedMaintainer.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/units/units.hpp>

#ifdef DYNAMO_DEBUG 
#include <boost/math/special_functions/fpclassify.hpp>
#endif

namespace dynamo {
  SysSchedMaintainer::SysSchedMaintainer(dynamo::Simulation* nSim, double ndt, std::string nName):
    System(nSim),
    periodt(ndt * nSim->units.unitTime())
  {
    dt = ndt * Sim->units.unitTime();
    sysName = nName;

    dout << "Periodic scheduler rebuild set for dt=" 
	 << ndt << std::endl;
  }

  void 
  SysSchedMaintainer::runEvent() const
  {
    double locdt = dt;
  
#ifdef DYNAMO_DEBUG 
    if (boost::math::isnan(dt))
      M_throw() << "A NAN system event time has been found";
#endif
    
    Sim->dSysTime += locdt;
    
    Sim->ptrScheduler->stream(locdt);
  
    //dynamics must be updated first
    Sim->stream(locdt);
  
    Sim->freestreamAcc += locdt;
  
    dt = periodt;
  
    Sim->ptrScheduler->rebuildList();
  }

  void 
  SysSchedMaintainer::initialise(size_t nID)
  { ID = nID; }

  void 
  SysSchedMaintainer::setdt(double ndt)
  { dt = ndt * Sim->units.unitTime(); }

  void 
  SysSchedMaintainer::increasedt(double ndt)
  { 
    dt += ndt * Sim->units.unitTime(); 
  }
}
