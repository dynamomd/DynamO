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

#include <dynamo/dynamics/systems/rescale.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/dynamics/units/units.hpp>
#include <dynamo/dynamics/BC/BC.hpp>
#include <dynamo/simulation/particle.hpp>
#include <dynamo/dynamics/species/species.hpp>
#include <dynamo/dynamics/NparticleEventData.hpp>
#include <dynamo/dynamics/ranges/include.hpp>
#include <dynamo/dynamics/liouvillean/liouvillean.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <boost/foreach.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <fstream>

namespace dynamo {
  SysRescale::SysRescale(const magnet::xml::Node& XML, dynamo::SimData* tmp): 
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

  SysRescale::SysRescale(dynamo::SimData* tmp, size_t frequency, std::string name, double kT):
    System(tmp),
    _frequency(frequency),
    _kT(kT),
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
  SysRescale::runEvent() const
  {
    double locdt = dt;

    Sim->dSysTime += locdt;
  
    Sim->ptrScheduler->stream(locdt);
  
    //dynamics must be updated first
    Sim->stream(locdt);
  
    ++Sim->eventCount;
    
    double currentkT(Sim->liouvillean->getkT()
		     / Sim->dynamics.units().unitEnergy());

    dout << "Rescaling kT " << currentkT 
	 << " To " << _kT / Sim->dynamics.units().unitEnergy() <<  std::endl;

    NEventData SDat;

    BOOST_FOREACH(const shared_ptr<Species>& species, Sim->species)
      BOOST_FOREACH(const unsigned long& partID, *species->getRange())
      SDat.L1partChanges.push_back(ParticleEventData(Sim->particleList[partID], *species, RESCALE));

    Sim->liouvillean->updateAllParticles();
    Sim->liouvillean->rescaleSystemKineticEnergy(_kT / currentkT);

    RealTime += (Sim->dSysTime - LastTime) / std::exp(0.5 * scaleFactor);

    LastTime = Sim->dSysTime;

    scaleFactor += std::log(currentkT);

    Sim->signalParticleUpdate(SDat);
  
    //Only 1ParticleEvents occur
    BOOST_FOREACH(const ParticleEventData& PDat, SDat.L1partChanges)
      Sim->ptrScheduler->fullUpdate(Sim->particleList[PDat.getParticle().getID()]);
  
    locdt += Sim->freestreamAcc;

    Sim->freestreamAcc = 0;

    BOOST_FOREACH(shared_ptr<OutputPlugin>& Ptr, Sim->outputPlugins)
      Ptr->eventUpdate(*this, SDat, locdt); 

    BOOST_FOREACH(shared_ptr<OutputPlugin>& Ptr, Sim->outputPlugins)
      Ptr->temperatureRescale(1.0/currentkT);

    dt = _timestep;
  
    Sim->ptrScheduler->rebuildList();
  }

  void 
  SysRescale::initialise(size_t nID)
  {
    ID = nID;

    dt = _timestep;

    if (_frequency != std::numeric_limits<size_t>::max())
      Sim->registerParticleUpdateFunc
	(magnet::function::MakeDelegate(this, &SysRescale::checker));
  
    dout << "Velocity rescaler initialising" << std::endl;
  }

  void 
  SysRescale::operator<<(const magnet::xml::Node& XML)
  {
    if (strcmp(XML.getAttribute("Type"),"Rescale"))
      M_throw() << "Attempting to load Rescale from " 
		<< XML.getAttribute("Type") << " entry"; 
  
    try {
      if (XML.hasAttribute("Freq"))
	_frequency = XML.getAttribute("Freq").as<size_t>();
    
      if (XML.hasAttribute("kT"))
	_kT = XML.getAttribute("kT").as<double>();
    
      _kT *= Sim->dynamics.units().unitEnergy();

      if (XML.hasAttribute("TimeStep"))
	_timestep = XML.getAttribute("TimeStep").as<double>();

      _timestep *= Sim->dynamics.units().unitTime();

      sysName = XML.getAttribute("Name");
    }
    catch (boost::bad_lexical_cast &)
      {
	M_throw() << "Failed a lexical cast in SysRescale";
      }
  }

  void 
  SysRescale::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::tag("System")
	<< magnet::xml::attr("Type") << "Rescale"
	<< magnet::xml::attr("kT") << _kT / Sim->dynamics.units().unitEnergy()
	<< magnet::xml::attr("Name") << sysName;
  
    if (_frequency != std::numeric_limits<size_t>::max())
      XML << magnet::xml::attr("Freq") << _frequency;

    if (_timestep != HUGE_VAL)
      XML << magnet::xml::attr("TimeStep") << _timestep / Sim->dynamics.units().unitTime();

    XML<< magnet::xml::endtag("System");
  }
}
