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

#include "rescale.hpp"
#include "../dynamics.hpp"
#include "../units/units.hpp"
#include "../BC/BC.hpp"
#include "../../simulation/particle.hpp"
#include "../species/species.hpp"
#include "../NparticleEventData.hpp"
#include "../ranges/include.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../../schedulers/scheduler.hpp"
#include <boost/foreach.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <fstream>

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
  Sim->dynamics.stream(locdt);
  
  ++Sim->eventCount;
  
  /////////Now the actual updates
  dout << "WARNING Rescaling kT to 1" << std::endl;
  
  double currentkT(Sim->dynamics.getLiouvillean().getkT()
		   / Sim->dynamics.units().unitEnergy());

  dout << "Current kT " << currentkT << std::endl;

  NEventData SDat;

  BOOST_FOREACH(const magnet::ClonePtr<Species>& species, Sim->dynamics.getSpecies())
    BOOST_FOREACH(const unsigned long& partID, *species->getRange())
    SDat.L1partChanges.push_back(ParticleEventData(Sim->particleList[partID], *species, RESCALE));

  Sim->dynamics.getLiouvillean().updateAllParticles();
  Sim->dynamics.getLiouvillean().rescaleSystemKineticEnergy(1.0/currentkT);

  RealTime += (Sim->dSysTime - LastTime) / std::exp(0.5 * scaleFactor);

  LastTime = Sim->dSysTime;

  scaleFactor += std::log(currentkT);

  Sim->signalParticleUpdate(SDat);
  
  //Only 1ParticleEvents occur
  BOOST_FOREACH(const ParticleEventData& PDat, SDat.L1partChanges)
    Sim->ptrScheduler->fullUpdate(PDat.getParticle());
  
  locdt += Sim->freestreamAcc;

  Sim->freestreamAcc = 0;

  BOOST_FOREACH(magnet::ClonePtr<OutputPlugin>& Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(*this, SDat, locdt); 

  BOOST_FOREACH(magnet::ClonePtr<OutputPlugin>& Ptr, Sim->outputPlugins)
    Ptr->temperatureRescale(1.0/currentkT);

  dt = _timestep;
  
  Sim->ptrScheduler->rebuildList();
}

void 
SysRescale::initialise(size_t nID)
{
  ID = nID;

  dt = HUGE_VAL;

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
SysRescale::outputXML(xml::XmlStream& XML) const
{
  XML << xml::tag("System")
      << xml::attr("Type") << "Rescale"
      << xml::attr("kT") << _kT / Sim->dynamics.units().unitEnergy()
      << xml::attr("Name") << sysName;
  
  if (_frequency != std::numeric_limits<size_t>::max())
    XML << xml::attr("Freq") << _frequency;

  if (_timestep != HUGE_VAL)
    XML << xml::attr("TimeStep") << _timestep / Sim->dynamics.units().unitTime();

  XML<< xml::endtag("System");
}
