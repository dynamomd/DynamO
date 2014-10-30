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
    _eventPeriod = 0;
    sysName = nName;

    dout << "Snapshot set for a period of " << _period / Sim->units.unitTime() << std::endl;
  }

  SysSnapshot::SysSnapshot(dynamo::Simulation* nSim, size_t nPeriod, std::string nName, std::string format, bool applyBC):
    System(nSim),
    _applyBC(applyBC),
    _format(format),
    _saveCounter(0)
  {
    _period = 0;
    dt = HUGE_VAL;
    _eventPeriod = nPeriod;
    sysName = nName;

    dout << "Snapshot set for a period of " << nPeriod << " events" << std::endl;
  }

  void
  SysSnapshot::eventCallback(const NEventData&)
  {
    if ((Sim->eventCount -_lastEventCount) >= _eventPeriod)
      {
	_lastEventCount = Sim->eventCount;
	dt = -HUGE_VAL;
	Sim->ptrScheduler->rebuildSystemEvents();
      }
  }
  
  NEventData
  SysSnapshot::runEvent()
  {
    if (_eventPeriod) 
      dt = HUGE_VAL;
    else
      dt += _period;

    Sim->dynamics->updateAllParticles();
  
    std::string filename = magnet::string::search_replace("Snapshot."+_format+".xml.bz2", "%COUNT", boost::lexical_cast<std::string>(_saveCounter));
    filename = magnet::string::search_replace(filename, "%ID", boost::lexical_cast<std::string>(Sim->simID));
    Sim->writeXMLfile(filename, _applyBC);
    
    dout << "Printing SNAPSHOT" << std::endl;
    
    filename = magnet::string::search_replace("Snapshot.output."+_format+".xml.bz2", "%COUNT", boost::lexical_cast<std::string>(_saveCounter++));
    filename = magnet::string::search_replace(filename, "%ID", boost::lexical_cast<std::string>(Sim->simID));
    Sim->outputData(filename);
    return NEventData();
  }

  void 
  SysSnapshot::initialise(size_t nID)
  { 
    ID = nID;
    _lastEventCount = Sim->eventCount;
    if (_eventPeriod)
      Sim->_sigParticleUpdate.connect<SysSnapshot, &SysSnapshot::eventCallback>(this);
  }
  
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
