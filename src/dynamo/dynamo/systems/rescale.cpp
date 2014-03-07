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

#include <dynamo/systems/rescale.hpp>

#include <dynamo/units/units.hpp>
#include <dynamo/BC/BC.hpp>
#include <dynamo/particle.hpp>
#include <dynamo/species/species.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/ranges/include.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <fstream>

namespace dynamo {
  SysRescale::SysRescale(const magnet::xml::Node& XML, dynamo::Simulation* tmp): 
    System(tmp),
    _frequency(std::numeric_limits<size_t>::max()),
    _kT(1),
    _timestep(HUGE_VAL),
    scaleFactor(1),
    LastTime(0),
    RealTime(0)
  {
    operator<<(XML);
    type = RESCALE;

    dout << "Velocity Rescaler Loaded" << std::endl;
  }

  SysRescale::SysRescale(dynamo::Simulation* tmp, size_t frequency, std::string name, double kT):
    System(tmp),
    _frequency(frequency),
    _kT(kT),
    _timestep(HUGE_VAL),
    scaleFactor(0),
    LastTime(0),
    RealTime(0)
  {
    type = RESCALE;
    sysName = name;

    dout << "Velocity Rescaler Loaded" << std::endl;
  }

  void 
  SysRescale::checker(const NEventData&)
  {
    if (!(Sim->eventCount % _frequency))
      {
	dt = 0;
	Sim->ptrScheduler->rebuildSystemEvents();
      }
  }

  void 
  SysRescale::runEvent()
  {
    double locdt = dt;

    Sim->systemTime += locdt;
  
    Sim->ptrScheduler->stream(locdt);
  
    //dynamics must be updated first
    Sim->stream(locdt);
  
    ++Sim->eventCount;
    
    double currentkT(Sim->dynamics->getkT()
		     / Sim->units.unitEnergy());

    dout << "Rescaling kT " << currentkT 
	 << " To " << _kT / Sim->units.unitEnergy() <<  std::endl;

    NEventData SDat;
    for (const shared_ptr<Species>& species : Sim->species)
      for (const unsigned long& partID : *species->getRange())
	SDat.L1partChanges.push_back(ParticleEventData(Sim->particles[partID], *species, RESCALE));
    
    Sim->dynamics->updateAllParticles();
    Sim->dynamics->rescaleSystemKineticEnergy(_kT / currentkT);
    
    //We must set the centre of mass velocity back to zero (assuming
    //this is the target velocity), otherwise it will drift with the
    //rescaling process
    Sim->setCOMVelocity();

    RealTime += (Sim->systemTime - LastTime) / std::exp(0.5 * scaleFactor);
    
    LastTime = Sim->systemTime;
    
    scaleFactor += std::log(currentkT);

    Sim->_sigParticleUpdate(SDat);
  
    //Only 1ParticleEvents occur
    for (const ParticleEventData& PDat : SDat.L1partChanges)
      Sim->ptrScheduler->fullUpdate(Sim->particles[PDat.getParticleID()]);
  
    for (shared_ptr<OutputPlugin>& Ptr : Sim->outputPlugins)
      Ptr->eventUpdate(*this, SDat, locdt); 

    dt = _timestep;

    Sim->ptrScheduler->rebuildList();
  }

  void 
  SysRescale::initialise(size_t nID)
  {
    ID = nID;

    dt = _timestep;

    if (_frequency != std::numeric_limits<size_t>::max())
      Sim->_sigParticleUpdate.connect<SysRescale, &SysRescale::checker>(this);
  
    dout << "Velocity rescaler initialising" << std::endl;
  }

  void 
  SysRescale::operator<<(const magnet::xml::Node& XML)
  {
    if (XML.hasAttribute("Freq"))
      _frequency = XML.getAttribute("Freq").as<size_t>();
    if (XML.hasAttribute("kT"))
      _kT = XML.getAttribute("kT").as<double>();
    _kT *= Sim->units.unitEnergy();
    if (XML.hasAttribute("TimeStep"))
      _timestep = XML.getAttribute("TimeStep").as<double>();
    _timestep *= Sim->units.unitTime();

    sysName = XML.getAttribute("Name");
  }

  void 
  SysRescale::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::tag("System")
	<< magnet::xml::attr("Type") << "Rescale"
	<< magnet::xml::attr("kT") << _kT / Sim->units.unitEnergy()
	<< magnet::xml::attr("Name") << sysName;
  
    if (_frequency != std::numeric_limits<size_t>::max())
      XML << magnet::xml::attr("Freq") << _frequency;

    if (_timestep != HUGE_VAL)
      XML << magnet::xml::attr("TimeStep") << _timestep / Sim->units.unitTime();

    XML<< magnet::xml::endtag("System");
  }
}
