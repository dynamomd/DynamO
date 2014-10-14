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

#include <dynamo/systems/francesco.hpp>
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

  SysFrancesco::SysFrancesco(const magnet::xml::Node& XML, dynamo::Simulation* tmp): 
    System(tmp),
    meanFreeTime(100000),
    Temp(Sim->units.unitEnergy()),
    sqrtTemp(std::sqrt(Sim->units.unitEnergy())),
    dimensions(NDIM),
    eventCount(0),
    lastlNColl(0)
  {
    dt = HUGE_VAL;
    operator<<(XML);
    type = GAUSSIAN;
  }

  SysFrancesco::SysFrancesco(dynamo::Simulation* nSim, double mft, double t, std::string nName):
    System(nSim),
    meanFreeTime(mft / nSim->N()),
    Temp(t),
    dimensions(NDIM),
    eventCount(0),
    lastlNColl(0),
    range(new IDRangeAll(Sim))
  {
    sysName = nName;
    type = GAUSSIAN;
  }

  void 
  SysFrancesco::runEvent()
  {
    Event event = getEvent();
    ++Sim->eventCount;
    ++eventCount;

#ifdef DYNAMO_DEBUG 
    if (std::isnan(event._dt))
      M_throw() << "A NAN system event time has been found";
#endif
    
    Sim->systemTime += event._dt;
    
    Sim->ptrScheduler->stream(event._dt);
  
    Sim->stream(event._dt);

    dt = getGhostt();

    size_t step = std::uniform_int_distribution<size_t>(0, range->size() - 1)(Sim->ranGenerator);
    Particle& part(Sim->particles[*(range->begin()+step)]);

    Sim->dynamics->updateParticle(part);
    ParticleEventData eventdata(part, *Sim->species(part), GAUSSIAN);

    //Locate surrounding particles, and calculate the average direction
    size_t n = 0;
    Vector avgV{0,0,0};
    std::unique_ptr<IDRange> ids(Sim->ptrScheduler->getParticleNeighbours(part));
    for (size_t ID2 : *ids)
      {
	auto& p2 = Sim->particles[ID2];
	Vector rij = part.getPosition() - p2.getPosition();
	Sim->BCs->applyBC(rij);
	if (rij.nrm2() > _R * _R) continue;
	Sim->dynamics->updateParticle(p2);
	avgV += p2.getVelocity().normal();
	++n;
      }
    avgV /= n;
    
    const double mass = Sim->species[eventdata.getSpeciesID()]->getMass(part);
    const double factor = sqrtTemp / std::sqrt(mass);
    //Get the new velocity
    std::normal_distribution<> norm_dist;
    double vel = std::abs(norm_dist(Sim->ranGenerator)) * factor;
    part.getVelocity() = avgV * vel;
    Sim->_sigParticleUpdate(eventdata);

    Sim->ptrScheduler->fullUpdate(part);
  
    for (shared_ptr<OutputPlugin>& Ptr : Sim->outputPlugins)
      Ptr->eventUpdate(event, eventdata);
  }

  void 
  SysFrancesco::initialise(size_t nID)
  {
    ID = nID;
    dt = getGhostt();
    sqrtTemp = sqrt(Temp);
    eventCount = 0;
    lastlNColl = 0;
    
    if (_R > Sim->ptrScheduler->getNeighbourhoodDistance())
      M_throw() << "The neighbourhood is too small for the R set in the Francesco System.";
  }

  void 
  SysFrancesco::operator<<(const magnet::xml::Node& XML)
  {
    meanFreeTime = XML.getAttribute("MFT").as<double>() * Sim->units.unitTime() / Sim->N();
    Temp = XML.getAttribute("Temperature").as<double>() * Sim->units.unitEnergy();
    sysName = XML.getAttribute("Name");

    _R = XML.getAttribute("R").as<double>() * Sim->units.unitLength();

    if (XML.hasAttribute("Dimensions"))
      dimensions = XML.getAttribute("Dimensions").as<size_t>();

    range = shared_ptr<IDRange>(IDRange::getClass(XML.getNode("IDRange"),Sim));
  }

  void 
  SysFrancesco::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::tag("System")
	<< magnet::xml::attr("Type") << "Francesco"
	<< magnet::xml::attr("Name") << sysName
	<< magnet::xml::attr("MFT") << meanFreeTime * Sim->N() / Sim->units.unitTime()
	<< magnet::xml::attr("Temperature") << Temp 
	<< magnet::xml::attr("R") << _R / Sim->units.unitLength() 
      / Sim->units.unitEnergy();
  
    if (dimensions != NDIM)
      XML << magnet::xml::attr("Dimensions") << dimensions;

    XML << range
	<< magnet::xml::endtag("System");
  }

  double 
  SysFrancesco::getGhostt() const
  { 
    return  - meanFreeTime * std::log(1.0 - std::uniform_real_distribution<>()(Sim->ranGenerator));
  }

  double 
  SysFrancesco::getReducedTemperature() const
  {
    return Temp / Sim->units.unitEnergy();
  }


  void 
  SysFrancesco::setReducedTemperature(double nT)
  {
    Temp = nT * Sim->units.unitEnergy(); 
  }
}
