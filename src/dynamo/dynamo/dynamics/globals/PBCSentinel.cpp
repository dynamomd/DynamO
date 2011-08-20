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

#include <dynamo/dynamics/globals/PBCSentinel.hpp>
#include <dynamo/dynamics/globals/globEvent.hpp>
#include <dynamo/dynamics/NparticleEventData.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <dynamo/dynamics/liouvillean/liouvillean.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <magnet/xmlreader.hpp>

#ifdef DYNAMO_DEBUG 
#include <boost/math/special_functions/fpclassify.hpp>
#endif

namespace dynamo {
  GPBCSentinel::GPBCSentinel(dynamo::SimData* nSim, const std::string& name):
    Global(nSim, "PBCSentinel"),
    maxintdist(0)
  {
    globName = name;
    dout << "PBCSentinel Loaded" << std::endl;
  }

  GPBCSentinel::GPBCSentinel(const magnet::xml::Node& XML, dynamo::SimData* ptrSim):
    Global(ptrSim, "PBCSentinel"),
    maxintdist(0)
  {
    operator<<(XML);

    dout << "PBCSentinel Loaded" << std::endl;
  }

  void 
  GPBCSentinel::initialise(size_t nID)
  {
    ID=nID;
  
    maxintdist = Sim->dynamics.getLongestInteraction();
  }

  void 
  GPBCSentinel::operator<<(const magnet::xml::Node& XML)
  {
    globName = XML.getAttribute("Name");	
  }

  GlobalEvent 
  GPBCSentinel::getEvent(const Particle& part) const
  {
    return GlobalEvent(part, Sim->dynamics.getLiouvillean().getPBCSentinelTime(part, maxintdist), 
		       VIRTUAL, *this);
  }

  void 
  GPBCSentinel::runEvent(const Particle& part, const double dt) const
  {
    GlobalEvent iEvent(part, dt, VIRTUAL, *this);

#ifdef DYNAMO_DEBUG 
    if (boost::math::isnan(iEvent.getdt()))
      M_throw() << "A NAN Interaction collision time has been found"
		<< iEvent.stringData(Sim);
  
    if (iEvent.getdt() == HUGE_VAL)
      M_throw() << "An infinite Interaction (not marked as NONE) collision time has been found\n"
		<< iEvent.stringData(Sim);
#endif

    Sim->dSysTime += iEvent.getdt();
    
    Sim->ptrScheduler->stream(iEvent.getdt());
  
    Sim->dynamics.stream(iEvent.getdt());

    Sim->dynamics.getLiouvillean().updateParticle(part);

#ifdef DYNAMO_DEBUG
    iEvent.addTime(Sim->freestreamAcc);
  
    Sim->freestreamAcc = 0;

    NEventData EDat(ParticleEventData(part, Sim->dynamics.getSpecies(part), VIRTUAL));

    Sim->signalParticleUpdate(EDat);

    BOOST_FOREACH(magnet::ClonePtr<OutputPlugin> & Ptr, Sim->outputPlugins)
      Ptr->eventUpdate(iEvent, EDat);
#else
    Sim->freestreamAcc += iEvent.getdt();
#endif

    Sim->ptrScheduler->fullUpdate(part);
  }

  void 
  GPBCSentinel::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::tag("Global") 
	<< magnet::xml::attr("Type") << "PBCSentinel"
	<< magnet::xml::attr("Name") << globName
	<< magnet::xml::endtag("Global");
  }
}
