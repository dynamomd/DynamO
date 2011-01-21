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

#include "schedMaintainer.hpp"
#include "../../base/is_simdata.hpp"
#include "../../schedulers/scheduler.hpp"
#include "../NparticleEventData.hpp"
#include "../units/units.hpp"

CSSchedMaintainer::CSSchedMaintainer(DYNAMO::SimData* nSim, double ndt, std::string nName):
  System(nSim),
  periodt(ndt * nSim->dynamics.units().unitTime())
{
  dt = ndt * Sim->dynamics.units().unitTime();
  sysName = nName;

  I_cout() << "Periodic scheduler rebuild set for dt=" 
	   << ndt;
}

void 
CSSchedMaintainer::runEvent() const
{
  double locdt = dt;
  
#ifdef DYNAMO_DEBUG 
  if (boost::math::isnan(dt))
    M_throw() << "A NAN system event time has been found";
#endif
    
  Sim->dSysTime += locdt;
    
  Sim->ptrScheduler->stream(locdt);
  
  //dynamics must be updated first
  Sim->dynamics.stream(locdt);
  
  Sim->freestreamAcc += locdt;
  
  dt = periodt;
  
  Sim->ptrScheduler->rebuildList();
}

void 
CSSchedMaintainer::initialise(size_t nID)
{ ID = nID; }

void 
CSSchedMaintainer::setdt(double ndt)
{ dt = ndt * Sim->dynamics.units().unitTime(); }

void 
CSSchedMaintainer::increasedt(double ndt)
{ 
  dt += ndt * Sim->dynamics.units().unitTime(); 
}
