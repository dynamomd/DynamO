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

#include <dynamo/systems/snapshot.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/outputplugins/tickerproperty/ticker.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <magnet/string/searchreplace.hpp>

namespace dynamo {
  SysSnapshot::SysSnapshot(dynamo::Simulation* nSim, double nPeriod, std::string nName, std::string format, bool applyBC):
    System(nSim),
    _applyBC(applyBC),
    _format(format),
    _saveCounter(0)
  {
    if (nPeriod <= 0.0)
      nPeriod = 1.0;

    nPeriod *= Sim->units.unitTime();

    dt = nPeriod;
    _period = nPeriod;

    sysName = nName;

    dout << "Snapshot set for a peroid of " 
	 << _period / Sim->units.unitTime() << std::endl;
  }

  void
  SysSnapshot::runEvent() const
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
  
    dt += _period;
  
    //This is done here as most ticker properties require it
    Sim->dynamics->updateAllParticles();

    for (shared_ptr<OutputPlugin>& Ptr : Sim->outputPlugins)
      Ptr->eventUpdate(*this, NEventData(), locdt);
  
    std::string filename = magnet::string::search_replace("Snapshot."+_format+".xml.bz2", "%COUNT", boost::lexical_cast<std::string>(_saveCounter));
    filename = magnet::string::search_replace(filename, "%ID", boost::lexical_cast<std::string>(Sim->simID));
    Sim->writeXMLfile(filename, _applyBC);
    
    dout << "Printing SNAPSHOT" << std::endl;
    
    filename = magnet::string::search_replace("Snapshot.output."+_format+".xml.bz2", "%COUNT", boost::lexical_cast<std::string>(_saveCounter++));
    filename = magnet::string::search_replace(filename, "%ID", boost::lexical_cast<std::string>(Sim->simID));
    Sim->outputData(filename);
  }

  void 
  SysSnapshot::initialise(size_t nID)
  { ID = nID; }

  void 
  SysSnapshot::setdt(double ndt)
  { 
    dt = ndt * Sim->units.unitTime(); 
  }

  void 
  SysSnapshot::increasedt(double ndt)
  { 
    dt += ndt * Sim->units.unitTime(); 
  }

  void 
  SysSnapshot::setTickerPeriod(const double& nP)
  { 
    dout << "Setting system ticker period to " 
	 << nP / Sim->units.unitTime() << std::endl;

    _period = nP; 

    dt = nP;

    if ((Sim->status >= INITIALISED) && Sim->endEventCount)
      Sim->ptrScheduler->rebuildSystemEvents();
  }
}
