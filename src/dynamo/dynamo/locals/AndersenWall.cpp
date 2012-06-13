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

#include <dynamo/locals/AndersenWall.hpp>
#include <dynamo/locals/localEvent.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/particle.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <magnet/overlap/cube_plane.hpp>
#include <magnet/xmlwriter.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_real.hpp>
#include <cmath>

namespace dynamo {
  LAndersenWall::LAndersenWall(const magnet::xml::Node& XML, dynamo::Simulation* ptrSim):
    Local(ptrSim, "GlobalAndersenWall"),
    sqrtT(1.0)
  {
    operator<<(XML);
  }

  LAndersenWall::LAndersenWall(dynamo::Simulation* nSim, double nsqrtT,
				 Vector  nnorm, Vector norigin, 
				 std::string nname, Range* nRange):
    Local(nRange, nSim, "AndersenWall"),
    vNorm(nnorm),
    vPosition(norigin),
    sqrtT(nsqrtT)
  {
    localName = nname;
  }

  LocalEvent 
  LAndersenWall::getEvent(const Particle& part) const
  {
#ifdef ISSS_DEBUG
    if (!Sim->dynamics->isUpToDate(part))
      M_throw() << "Particle is not up to date";
#endif

    return LocalEvent(part, Sim->dynamics->getWallCollision(part, vPosition, vNorm), WALL, *this);
  }

  void
  LAndersenWall::runEvent(Particle& part, const LocalEvent& iEvent) const
  {
    ++Sim->eventCount;
  
    NEventData EDat(Sim->dynamics->runAndersenWallCollision
		    (part, vNorm, sqrtT));
  
    Sim->signalParticleUpdate(EDat);
  
    Sim->ptrScheduler->fullUpdate(part);
  
    BOOST_FOREACH(shared_ptr<OutputPlugin> & Ptr, Sim->outputPlugins)
      Ptr->eventUpdate(iEvent, EDat);
  }

  bool 
  LAndersenWall::isInCell(const Vector & Origin, 
			   const Vector & CellDim) const
  {
    return magnet::overlap::cube_plane(Origin, CellDim, vPosition, vNorm);
  }

  void 
  LAndersenWall::initialise(size_t nID)
  {
    ID = nID;
  }

  void 
  LAndersenWall::operator<<(const magnet::xml::Node& XML)
  {
    range = shared_ptr<Range>(Range::getClass(XML,Sim));
  
    try {
    
      sqrtT = sqrt(XML.getAttribute("Temperature").as<double>() 
		   * Sim->units.unitEnergy());

      localName = XML.getAttribute("Name");
      vNorm << XML.getNode("Norm");
      vNorm /= vNorm.nrm();
      vPosition << XML.getNode("Origin");
      vPosition *= Sim->units.unitLength();

    } 
    catch (boost::bad_lexical_cast &)
      {
	M_throw() << "Failed a lexical cast in LAndersenWall";
      }
  }

  void
  LAndersenWall::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "AndersenWall"
	<< magnet::xml::attr("Name") << localName
	<< magnet::xml::attr("Temperature") << sqrtT * sqrtT 
      / Sim->units.unitEnergy()
	<< range
	<< magnet::xml::tag("Norm")
	<< vNorm
	<< magnet::xml::endtag("Norm")
	<< magnet::xml::tag("Origin")
	<< vPosition / Sim->units.unitLength()
	<< magnet::xml::endtag("Origin");
  }
}
