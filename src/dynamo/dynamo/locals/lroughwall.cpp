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

#include <dynamo/locals/lroughwall.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/BC/BC.hpp>
#include <dynamo/locals/localEvent.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <magnet/overlap/cube_plane.hpp>

namespace dynamo {
  LRoughWall::LRoughWall(dynamo::Simulation* nSim, double ne, double net, double nr, Vector  nnorm, 
			 Vector  norigin, std::string nname, 
			 IDRange* nRange, bool nrender):
    Local(nRange, nSim, "LocalRoughWall"),
    vNorm(nnorm),
    vPosition(norigin),
    e(ne),
    et(net),
    r(nr),
    render(nrender)
  {
    localName = nname;
  }

  LRoughWall::LRoughWall(const magnet::xml::Node& XML, dynamo::Simulation* tmp):
    Local(tmp, "LocalRoughWall")
  {
    operator<<(XML);
  }

  LocalEvent 
  LRoughWall::getEvent(const Particle& part) const
  {
#ifdef ISSS_DEBUG
    if (!Sim->dynamics->isUpToDate(part))
      M_throw() << "Particle is not up to date";
#endif

    return LocalEvent(part, Sim->dynamics->getPlaneEvent(part, vPosition, vNorm, r), WALL, *this);
  }

  void
  LRoughWall::runEvent(Particle& part, const LocalEvent& iEvent) const
  {
    ++Sim->eventCount;

    //Run the collision and catch the data
    NEventData EDat(Sim->dynamics->runRoughWallCollision
		    (part, vNorm, e, et, r));

    (*Sim->_sigParticleUpdate)(EDat);

    //Now we're past the event update the scheduler and plugins
    Sim->ptrScheduler->fullUpdate(part);
  
    for (shared_ptr<OutputPlugin> & Ptr : Sim->outputPlugins)
      Ptr->eventUpdate(iEvent, EDat);
  }

  void 
  LRoughWall::operator<<(const magnet::xml::Node& XML)
  {
    range = shared_ptr<IDRange>(IDRange::getClass(XML.getNode("IDRange"),Sim));
    e = XML.getAttribute("Elasticity").as<double>();
    et = XML.getAttribute("TangentialElasticity").as<double>();
    r = XML.getAttribute("Radius").as<double>() * Sim->units.unitLength();
    render = XML.getAttribute("Render").as<double>();
    localName = XML.getAttribute("Name");
    vNorm << XML.getNode("Norm");
    vNorm /= vNorm.nrm();
    vPosition << XML.getNode("Origin");
    vPosition *= Sim->units.unitLength();
  }

  void 
  LRoughWall::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "RoughWall" 
	<< magnet::xml::attr("Name") << localName
	<< magnet::xml::attr("Elasticity") << e
	<< magnet::xml::attr("TangentialElasticity") << et
	<< magnet::xml::attr("Radius") << r / Sim->units.unitLength()
	<< magnet::xml::attr("Render") << render
	<< range
	<< magnet::xml::tag("Norm")
	<< vNorm
	<< magnet::xml::endtag("Norm")
	<< magnet::xml::tag("Origin")
	<< vPosition / Sim->units.unitLength()
	<< magnet::xml::endtag("Origin");
  }

  bool 
  LRoughWall::validateState(const Particle& part, bool textoutput) const
  {
    Vector pos(part.getPosition() - vPosition);
    Sim->BCs->applyBC(pos);
    
    double overlap = r - std::abs(pos | vNorm);
    
    if (overlap > 0)
      {
	if (textoutput)
	  derr << "Particle " << part.getID() << " is " << overlap / Sim->units.unitLength() << " far into the wall."
	       << "\nWall Pos = " << Vector(vPosition / Sim->units.unitLength()).toString()
	       << ", Normal = " << vNorm.toString() << ", r = " << r / Sim->units.unitLength()
	       << std::endl;
	return true;
      }
    return false;
  }
}

