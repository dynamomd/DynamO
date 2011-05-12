/*  dynamo:- Event driven molecular dynamics simulator 
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

#include "snapshot.hpp"
#include "../../base/is_simdata.hpp"
#include "../NparticleEventData.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../../outputplugins/tickerproperty/ticker.hpp"
#include "../units/units.hpp"
#include "../../schedulers/scheduler.hpp"

SSnapshot::SSnapshot(dynamo::SimData* nSim, double nPeriod, std::string nName):
  System(nSim),
  _applyBC(false),
  _saveCounter(0)

{
  if (nPeriod <= 0.0)
    nPeriod = 1.0;

  nPeriod *= Sim->dynamics.units().unitTime();

  dt = nPeriod;
  _period = nPeriod;

  sysName = nName;

  I_cout() << "Snapshot set for a peroid of " 
	   << _period / Sim->dynamics.units().unitTime();
}

void
SSnapshot::runEvent() const
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
  
  dt += _period;
  
  locdt += Sim->freestreamAcc;
  Sim->freestreamAcc = 0;

  //This is done here as most ticker properties require it
  Sim->dynamics.getLiouvillean().updateAllParticles();

  BOOST_FOREACH(magnet::ClonePtr<OutputPlugin>& Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(*this, NEventData(), locdt);
  
  std::string filename = dynamo::searchReplace("Snapshot.%i.xml.bz2", "%i", boost::lexical_cast<std::string>(_saveCounter++));
  Sim->writeXMLfile(filename, _applyBC);
}

void 
SSnapshot::initialise(size_t nID)
{ ID = nID; }

void 
SSnapshot::setdt(double ndt)
{ 
  dt = ndt * Sim->dynamics.units().unitTime(); 
}

void 
SSnapshot::increasedt(double ndt)
{ 
  dt += ndt * Sim->dynamics.units().unitTime(); 
}

void 
SSnapshot::setTickerPeriod(const double& nP)
{ 
  I_cout() << "Setting system ticker period to " 
	   << nP / Sim->dynamics.units().unitTime();

  _period = nP; 

  dt = nP;

  if ((Sim->status >= INITIALISED) && Sim->endEventCount)
    Sim->ptrScheduler->rebuildSystemEvents();
}
