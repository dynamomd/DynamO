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

#include <dynamo/globals/francesco.hpp>
#include <dynamo/globals/globEvent.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
  GFrancesco::GFrancesco(const magnet::xml::Node& XML, dynamo::Simulation* ptrSim):
    Global(ptrSim, "Francesco")
  {
    operator<<(XML);
  }

  void 
  GFrancesco::initialise(size_t nID)
  {
    Global::initialise(nID);
    _eventTimes.clear();
    _eventTimes.resize(Sim->particles.size(), HUGE_VAL);
    Sim->_sigParticleUpdate.connect<GFrancesco, &GFrancesco::particlesUpdated>(this);
  }

  void 
  GFrancesco::particlesUpdated(const NEventData& PDat)
  {
    for (const PairEventData& pdat : PDat.L2partChanges)
      {
	_eventTimes[pdat.particle1_.getParticleID()] = - _MFT * std::log(1.0 - std::uniform_real_distribution<>()(Sim->ranGenerator)) + Sim->systemTime;
	_eventTimes[pdat.particle2_.getParticleID()] = - _MFT * std::log(1.0 - std::uniform_real_distribution<>()(Sim->ranGenerator)) + Sim->systemTime;
      }
  }
  void 
  GFrancesco::outputXML(magnet::xml::XmlStream& XML) const {
    XML << magnet::xml::tag("Global")
	<< magnet::xml::attr("Type") << "Francesco"
	<< magnet::xml::attr("Name") << globName
	<< magnet::xml::attr("MFT") << _MFT / Sim->units.unitTime()
	<< magnet::xml::attr("Temperature") << _T / Sim->units.unitEnergy()
	<< range
    	<< magnet::xml::endtag("Global");
  }

  void 
  GFrancesco::operator<<(const magnet::xml::Node& XML)
  {
    globName = XML.getAttribute("Name");
    _MFT = XML.getAttribute("MFT").as<double>() * Sim->units.unitTime();
    _T = XML.getAttribute("Temperature").as<double>() * Sim->units.unitEnergy();
    range = shared_ptr<IDRange>(IDRange::getClass(XML.getNode("IDRange"), Sim));
  }

  GlobalEvent 
  GFrancesco::getEvent(const Particle& part) const {
    return GlobalEvent(part, _eventTimes[part] - Sim->systemTime, GAUSSIAN, *this);
  }

  void 
  GFrancesco::runEvent(Particle& part, const double)
  {
    const double dt = _eventTimes[part] - Sim->systemTime;
    _eventTimes[part] = HUGE_VAL;
    GlobalEvent iEvent(part, dt, GAUSSIAN, *this);
    
    Sim->systemTime += dt;
    Sim->ptrScheduler->stream(dt);
    Sim->stream(dt);

    Sim->dynamics->updateParticle(part);
    NEventData EDat(ParticleEventData(part, *Sim->species(part), GAUSSIAN));
    
    //Kill the rotational motion
    Sim->dynamics->getRotData(part).angularVelocity = Vector(0,0,0);
    //Reassign the linear motion
    //const double mass = Sim->species[part]->getMass(part);
    //std::normal_distribution<> norm_dist(0, std::sqrt(_T / mass));
    //double vel = std::abs(norm_dist(Sim->ranGenerator));
    part.getVelocity() = Sim->dynamics->getRotData(part).orientation * magnet::math::Quaternion::initialDirector();

    Sim->_sigParticleUpdate(EDat);
    for (shared_ptr<OutputPlugin> & Ptr : Sim->outputPlugins)
      Ptr->eventUpdate(iEvent, EDat);
    Sim->ptrScheduler->fullUpdate(part);
  }
}
