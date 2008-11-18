/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

CSSchedMaintainer::CSSchedMaintainer(DYNAMO::SimData* nSim, Iflt ndt, std::string nName):
  CSystem(nSim),
  periodt(ndt * nSim->Dynamics.units().unitTime())
{
  dt = ndt * Sim->Dynamics.units().unitTime();
  sysName = nName;

  I_cout() << "Periodic scheduler rebuild set for dt=" 
	   << ndt;
}

void 
CSSchedMaintainer::stream(Iflt ndt)
{
  dt -= ndt;
}

void 
CSSchedMaintainer::runEvent() const
{
  Iflt locdt = dt;
  
#ifdef DYNAMO_DEBUG 
  if (isnan(dt))
    D_throw() << "A NAN system event time has been found";
#endif
    
  Sim->dSysTime += locdt;
    
  Sim->ptrScheduler->stream(locdt);
  
  //dynamics must be updated first
  Sim->Dynamics.stream(locdt);
  
  dt = periodt;

  //Run the collision and catch the data
  
  BOOST_FOREACH(smrtPlugPtr<COutputPlugin>& Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(static_cast<const CSystem&>(*this), CNParticleData(), locdt);

  Sim->ptrScheduler->rebuildList();
}

void 
CSSchedMaintainer::initialise(size_t nID)
{ ID = nID; }

void 
CSSchedMaintainer::setdt(Iflt ndt)
{ dt = ndt * Sim->Dynamics.units().unitTime(); }

void 
CSSchedMaintainer::increasedt(Iflt ndt)
{ 
  dt += ndt * Sim->Dynamics.units().unitTime(); 
}
