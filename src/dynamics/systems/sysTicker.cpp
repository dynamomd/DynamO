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

#include "sysTicker.hpp"
#include "../../base/is_simdata.hpp"
#include "../NparticleEventData.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../../outputplugins/tickerproperty/ticker.hpp"
#include "../units/units.hpp"
#include "../../schedulers/scheduler.hpp"

CSTicker::CSTicker(DYNAMO::SimData* nSim, Iflt nPeriod, std::string nName):
  CSystem(nSim)
{
  if (nPeriod <= 0.0)
    nPeriod = Sim->Dynamics.units().unitTime();

  dt = nPeriod;
  period = nPeriod;

  sysName = nName;

  I_cout() << "System ticker set for a peroid of " 
	   << nPeriod / Sim->Dynamics.units().unitTime();
}

void
CSTicker::runEvent() const
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
  
  dt += period;
  
  locdt += Sim->freestreamAcc;
  Sim->freestreamAcc = 0;

  //This is done here as most ticker properties require it
  Sim->Dynamics.getLiouvillean().updateAllParticles();

  {
    COPTicker* ptr = NULL;
    BOOST_FOREACH(smrtPlugPtr<COutputPlugin>& Ptr, Sim->outputPlugins)
      {
	ptr = dynamic_cast<COPTicker*>(Ptr.get_ptr());
	if (ptr != NULL)
	  ptr->ticker();
      }
  }

  BOOST_FOREACH(smrtPlugPtr<COutputPlugin>& Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(*this, CNParticleData(), locdt);
}

void 
CSTicker::initialise(size_t nID)
{ ID = nID; }

void 
CSTicker::setdt(Iflt ndt)
{ 
  dt = ndt * Sim->Dynamics.units().unitTime(); 
}

void 
CSTicker::increasedt(Iflt ndt)
{ 
  dt += ndt * Sim->Dynamics.units().unitTime(); 
}

void 
CSTicker::setTickerPeriod(const Iflt& nP)
{ 
  I_cout() << "Setting system ticker period to " 
	   << nP / Sim->Dynamics.units().unitTime();

  period = nP; 

  dt = nP;

  if ((Sim->status >= INITIALISED) && Sim->lMaxNColl)
    Sim->ptrScheduler->rebuildSystemEvents();
}
