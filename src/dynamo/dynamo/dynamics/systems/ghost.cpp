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

#include <dynamo/dynamics/systems/ghost.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/dynamics/units/units.hpp>
#include <dynamo/dynamics/BC/BC.hpp>
#include <dynamo/simulation/particle.hpp>
#include <dynamo/dynamics/species/species.hpp>
#include <dynamo/dynamics/NparticleEventData.hpp>
#include <dynamo/dynamics/ranges/include.hpp>
#include <dynamo/dynamics/liouvillean/liouvillean.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random/uniform_int.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

#ifdef DYNAMO_DEBUG 
#include <boost/math/special_functions/fpclassify.hpp>
#endif

namespace dynamo {

  SysAndersen::SysAndersen(const magnet::xml::Node& XML, dynamo::SimData* tmp): 
    System(tmp),
    meanFreeTime(100000),
    Temp(Sim->dynamics.units().unitEnergy()),
    sqrtTemp(std::sqrt(Sim->dynamics.units().unitEnergy())),
    tune(false),
    setPoint(0.05),
    eventCount(0),
    lastlNColl(0),
    setFrequency(100)
  {
    dt = HUGE_VAL;
    operator<<(XML);
    type = GAUSSIAN;
  }

  SysAndersen::SysAndersen(dynamo::SimData* nSim, double mft, double t, 
		       std::string nName):
    System(nSim),
    meanFreeTime(mft),
    Temp(t),
    tune(true),
    setPoint(0.05),
    eventCount(0),
    lastlNColl(0),
    setFrequency(100),
    range(new RAll(Sim))
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
    if (boost::math::isnan(locdt))
      M_throw() << "A NAN system event time has been found";
#endif
    
    Sim->dSysTime += locdt;
    
    Sim->ptrScheduler->stream(locdt);
  
    Sim->dynamics.stream(locdt);

    locdt +=  Sim->freestreamAcc;
    Sim->freestreamAcc = 0;

    dt = getGhostt();

    unsigned int step = boost::variate_generator
      <dynamo::baseRNG&, boost::uniform_int<unsigned int> >
      (Sim->ranGenerator, 
       boost::uniform_int<unsigned int>(0, range->size() - 1))();

    const Particle& part(Sim->particleList[*(range->begin()+step)]);

    //Run the collision and catch the data
    NEventData SDat(Sim->dynamics.getLiouvillean().randomGaussianEvent
		    (part, sqrtTemp));
  
    Sim->signalParticleUpdate(SDat);

    Sim->ptrScheduler->fullUpdate(part);
  
    BOOST_FOREACH(shared_ptr<OutputPlugin>& Ptr, Sim->outputPlugins)
      Ptr->eventUpdate(*this, SDat, locdt);

  }

  void 
  SysAndersen::initialise(size_t nID)
  {
    ID = nID;
    meanFreeTime /= Sim->N;
    dt = getGhostt();
    sqrtTemp = sqrt(Temp);
  }

  void 
  SysAndersen::operator<<(const magnet::xml::Node& XML)
  {
    if (strcmp(XML.getAttribute("Type"),"Andersen"))
      M_throw() << "Attempting to load Andersen from non Andersen entry"; 
  
    try {
      meanFreeTime = XML.getAttribute("MFT").as<double>() * Sim->dynamics.units().unitTime();
      Temp = XML.getAttribute("Temperature").as<double>() * Sim->dynamics.units().unitEnergy();
      sysName = XML.getAttribute("Name");

      if (XML.hasAttribute("SetFrequency") && XML.hasAttribute("SetPoint"))
	{
	  tune = true;
	  setFrequency = XML.getAttribute("SetFrequency").as<unsigned long long>();
	  setPoint = boost::lexical_cast<double>(XML.getAttribute("SetPoint"));
	}

      range = shared_ptr<Range>(Range::getClass(XML,Sim));
    }
    catch (boost::bad_lexical_cast &)
      {
	M_throw() << "Failed a lexical cast in CGGlobal";
      }
  }

  void 
  SysAndersen::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::tag("System")
	<< magnet::xml::attr("Type") << "Andersen"
	<< magnet::xml::attr("Name") << sysName
	<< magnet::xml::attr("MFT") << meanFreeTime
      * Sim->N
      / Sim->dynamics.units().unitTime()
	<< magnet::xml::attr("Temperature") << Temp 
      / Sim->dynamics.units().unitEnergy();
  
    if (tune)
      XML << magnet::xml::attr("SetPoint") << setPoint
	  << magnet::xml::attr("SetFrequency") << setFrequency;
  
    XML << range
	<< magnet::xml::endtag("System");
  }

  double 
  SysAndersen::getGhostt() const
  { 
    return  - meanFreeTime * log(1 - Sim->uniform_sampler());
  }

  double 
  SysAndersen::getReducedTemperature() const
  {
    return Temp / Sim->dynamics.units().unitEnergy();
  }


  void 
  SysAndersen::setReducedTemperature(double nT)
  {
    Temp = nT * Sim->dynamics.units().unitEnergy(); 
  }
}
