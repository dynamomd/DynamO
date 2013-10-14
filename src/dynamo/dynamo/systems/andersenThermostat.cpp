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

#include <dynamo/systems/andersenThermostat.hpp>

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

namespace dynamo {

  SysAndersen::SysAndersen(const magnet::xml::Node& XML, dynamo::Simulation* tmp): 
    System(tmp),
    meanFreeTime(100000),
    Temp(Sim->units.unitEnergy()),
    sqrtTemp(std::sqrt(Sim->units.unitEnergy())),
    tune(false),
    dimensions(NDIM),
    setPoint(0.05),
    eventCount(0),
    lastlNColl(0),
    setFrequency(100)
  {
    dt = HUGE_VAL;
    operator<<(XML);
    type = GAUSSIAN;
  }

  SysAndersen::SysAndersen(dynamo::Simulation* nSim, double mft, double t, 
		       std::string nName):
    System(nSim),
    meanFreeTime(mft),
    Temp(t),
    tune(true),
    dimensions(NDIM),
    setPoint(0.05),
    eventCount(0),
    lastlNColl(0),
    setFrequency(100),
    range(new IDRangeAll(Sim))
  {
    sysName = nName;
    type = GAUSSIAN;
  }

  void 
  SysAndersen::runEvent() const
  {
    ++Sim->eventCount;
    ++eventCount;

    if (tune && (eventCount > setFrequency))
      {
	meanFreeTime *= static_cast<double>(eventCount)
	  / ((Sim->eventCount - lastlNColl) * setPoint);

	lastlNColl = Sim->eventCount;
	eventCount = 0;
      }

    double locdt = dt;
  
#ifdef DYNAMO_DEBUG 
    if (std::isnan(locdt))
      M_throw() << "A NAN system event time has been found";
#endif
    
    Sim->systemTime += locdt;
    
    Sim->ptrScheduler->stream(locdt);
  
    Sim->stream(locdt);

    dt = getGhostt();

    size_t step = std::uniform_int_distribution<size_t>(0, range->size() - 1)(Sim->ranGenerator);

    Particle& part(Sim->particles[*(range->begin()+step)]);

    //Run the collision and catch the data
    NEventData SDat(Sim->dynamics->randomGaussianEvent(part, sqrtTemp, dimensions));
  
    Sim->_sigParticleUpdate(SDat);

    Sim->ptrScheduler->fullUpdate(part);
  
    for (shared_ptr<OutputPlugin>& Ptr : Sim->outputPlugins)
      Ptr->eventUpdate(*this, SDat, locdt);
  }

  void 
  SysAndersen::initialise(size_t nID)
  {
    ID = nID;
    meanFreeTime /= Sim->N();
    dt = getGhostt();
    sqrtTemp = sqrt(Temp);
  }

  void 
  SysAndersen::operator<<(const magnet::xml::Node& XML)
  {
    meanFreeTime = XML.getAttribute("MFT").as<double>() * Sim->units.unitTime();
    Temp = XML.getAttribute("Temperature").as<double>() * Sim->units.unitEnergy();
    sysName = XML.getAttribute("Name");

    if (XML.hasAttribute("Dimensions"))
      dimensions = XML.getAttribute("Dimensions").as<size_t>();

    if (XML.hasAttribute("SetFrequency") && XML.hasAttribute("SetPoint"))
      {
	tune = true;
	setFrequency = XML.getAttribute("SetFrequency").as<size_t>();
	setPoint = XML.getAttribute("SetPoint").as<double>();
      }

    range = shared_ptr<IDRange>(IDRange::getClass(XML.getNode("IDRange"),Sim));
  }

  void 
  SysAndersen::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::tag("System")
	<< magnet::xml::attr("Type") << "Andersen"
	<< magnet::xml::attr("Name") << sysName
	<< magnet::xml::attr("MFT") << meanFreeTime
      * Sim->N()
      / Sim->units.unitTime()
	<< magnet::xml::attr("Temperature") << Temp 
      / Sim->units.unitEnergy();
  
    if (tune)
      XML << magnet::xml::attr("SetPoint") << setPoint
	  << magnet::xml::attr("SetFrequency") << setFrequency;
  
    if (dimensions != NDIM)
      XML << magnet::xml::attr("Dimensions") << dimensions;

    XML << range
	<< magnet::xml::endtag("System");
  }

  double 
  SysAndersen::getGhostt() const
  { 
    return  - meanFreeTime * std::log(1.0 - std::uniform_real_distribution<>()(Sim->ranGenerator));
  }

  double 
  SysAndersen::getReducedTemperature() const
  {
    return Temp / Sim->units.unitEnergy();
  }


  void 
  SysAndersen::setReducedTemperature(double nT)
  {
    Temp = nT * Sim->units.unitEnergy(); 
  }
}
