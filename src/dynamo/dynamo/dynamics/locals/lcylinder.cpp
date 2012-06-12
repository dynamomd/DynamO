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

#include <dynamo/dynamics/locals/lcylinder.hpp>
#include <dynamo/dynamics/liouvillean/liouvillean.hpp>
#include <dynamo/dynamics/locals/localEvent.hpp>
#include <dynamo/dynamics/NparticleEventData.hpp>
#include <dynamo/dynamics/units/units.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>

namespace dynamo {
  LCylinder::LCylinder(dynamo::SimData* nSim, double ne, Vector  nnorm, 
			 Vector  norigin, double nr, std::string nname, 
			 Range* nRange, bool nrender):
    Local(nRange, nSim, "CylinderWall"),
    vNorm(nnorm),
    vPosition(norigin),
    e(ne),
    radius(nr),
    render(nrender)
  {
    localName = nname;
  }

  LCylinder::LCylinder(const magnet::xml::Node& XML, dynamo::SimData* tmp):
    Local(tmp, "CylinderWall")
  {
    operator<<(XML);
  }

  LocalEvent 
  LCylinder::getEvent(const Particle& part) const
  {
#ifdef ISSS_DEBUG
    if (!Sim->liouvillean->isUpToDate(part))
      M_throw() << "Particle is not up to date";
#endif

    return LocalEvent(part, Sim->liouvillean->getCylinderWallCollision
		      (part, vPosition, vNorm, radius), WALL, *this);
  }

  void
  LCylinder::runEvent(Particle& part, const LocalEvent& iEvent) const
  {
    ++Sim->eventCount;

    //Run the collision and catch the data
    NEventData EDat(Sim->liouvillean->runCylinderWallCollision
		    (part, vPosition, vNorm, e));

    Sim->signalParticleUpdate(EDat);

    //Now we're past the event update the scheduler and plugins
    Sim->ptrScheduler->fullUpdate(part);
  
    BOOST_FOREACH(shared_ptr<OutputPlugin> & Ptr, Sim->outputPlugins)
      Ptr->eventUpdate(iEvent, EDat);
  }

  bool 
  LCylinder::isInCell(const Vector & Origin, const Vector& CellDim) const
  {
    return true;
  }

  void 
  LCylinder::initialise(size_t nID)
  {
    ID = nID;
  }

  void 
  LCylinder::operator<<(const magnet::xml::Node& XML)
  {
    range = shared_ptr<Range>(Range::getClass(XML,Sim));
  
    try {
      e = XML.getAttribute("Elasticity").as<double>();

      radius = XML.getAttribute("Radius").as<double>() * Sim->dynamics.units().unitLength();

      render = XML.getAttribute("Render").as<bool>();
      magnet::xml::Node xBrowseNode = XML.getNode("Norm");
      localName = XML.getAttribute("Name");
      vNorm << xBrowseNode;
      vNorm /= vNorm.nrm();
      xBrowseNode = XML.getNode("Origin");
      vPosition << xBrowseNode;
      vPosition *= Sim->dynamics.units().unitLength();
    } 
    catch (boost::bad_lexical_cast &)
      {
	M_throw() << "Failed a lexical cast in LCylinder";
      }
  }

  void 
  LCylinder::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "CylinderWall" 
	<< magnet::xml::attr("Name") << localName
	<< magnet::xml::attr("Elasticity") << e
	<< magnet::xml::attr("Radius") << radius / Sim->dynamics.units().unitLength()
	<< magnet::xml::attr("Render") << render
	<< range
	<< magnet::xml::tag("Norm")
	<< vNorm
	<< magnet::xml::endtag("Norm")
	<< magnet::xml::tag("Origin")
	<< vPosition / Sim->dynamics.units().unitLength()
	<< magnet::xml::endtag("Origin");
  }
}

