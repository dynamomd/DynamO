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

#include <dynamo/dynamics/locals/ldblwall.hpp>
#include <dynamo/dynamics/BC/BC.hpp>
#include <dynamo/dynamics/liouvillean/liouvillean.hpp>
#include <dynamo/dynamics/locals/localEvent.hpp>
#include <dynamo/dynamics/NparticleEventData.hpp>
#include <dynamo/dynamics/units/units.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <magnet/overlap/cube_plane.hpp>

namespace dynamo {
  CLDblWall::CLDblWall(dynamo::SimData* nSim, double ne, Vector  nnorm, 
		       Vector  norigin, std::string nname, CRange* nRange):
    Local(nRange, nSim, "LocalWall"),
    vNorm(nnorm),
    vPosition(norigin),
    e(ne)
  {
    localName = nname;
  }

  CLDblWall::CLDblWall(const magnet::xml::Node& XML, dynamo::SimData* tmp):
    Local(tmp, "LocalDoubleWall")
  {
    operator<<(XML);
  }

  LocalEvent 
  CLDblWall::getEvent(const Particle& part) const
  {
#ifdef ISSS_DEBUG
    if (!Sim->dynamics.getLiouvillean().isUpToDate(part))
      M_throw() << "Particle is not up to date";
#endif

    if (part.getID() == lastID) return LocalEvent(part, HUGE_VAL, NONE, *this);
  
    Vector rij = part.getPosition() - vPosition;
    Sim->dynamics.BCs().applyBC(rij);

    Vector norm(vNorm);
    if ((norm | rij) < 0)
      norm *= -1;

    return LocalEvent(part, Sim->dynamics.getLiouvillean().getWallCollision
		      (part, vPosition, norm), WALL, *this);
  }

  void
  CLDblWall::runEvent(const Particle& part, const LocalEvent& iEvent) const
  {
    ++Sim->eventCount;

    Vector norm = vNorm;

    Vector rij = part.getPosition() - vPosition;
    Sim->dynamics.BCs().applyBC(rij);
  
    if ((norm | rij) < 0)
      norm *= -1;

    //Run the collision and catch the data
    NEventData EDat(Sim->dynamics.getLiouvillean().runWallCollision
		    (part, norm, e));

    Sim->signalParticleUpdate(EDat);

    //Must do this after the signal is run
    lastID = part.getID();

    //Now we're past the event update the scheduler and plugins
    Sim->ptrScheduler->fullUpdate(part);
  
    BOOST_FOREACH(std::tr1::shared_ptr<OutputPlugin> & Ptr, Sim->outputPlugins)
      Ptr->eventUpdate(iEvent, EDat);
  }

  bool 
  CLDblWall::isInCell(const Vector & Origin, const Vector& CellDim) const
  {
    return magnet::overlap::cube_plane(Origin, CellDim, vPosition, vNorm);
  }

  void 
  CLDblWall::initialise(size_t nID)
  {
    ID = nID;
    lastID = std::numeric_limits<size_t>::max();

    Sim->registerParticleUpdateFunc
      (magnet::function::MakeDelegate(this, &CLDblWall::particleUpdate));

  }

  void
  CLDblWall::particleUpdate(const NEventData& PDat) const
  {
    if (lastID == std::numeric_limits<size_t>::max()) return;

    BOOST_FOREACH(const ParticleEventData& pdat, PDat.L1partChanges)
      if (pdat.getParticle().getID() == lastID)
	{
	  lastID = std::numeric_limits<size_t>::max();
	  return;
	}
  
    BOOST_FOREACH(const PairEventData& pdat, PDat.L2partChanges)
      if (pdat.particle1_.getParticle().getID() == lastID 
	  || pdat.particle2_.getParticle().getID() == lastID)
	{
	  lastID = std::numeric_limits<size_t>::max();
	  return;
	}
  }

  void 
  CLDblWall::operator<<(const magnet::xml::Node& XML)
  {
    range = std::tr1::shared_ptr<CRange>(CRange::getClass(XML,Sim));
  
    try {
      e = XML.getAttribute("Elasticity").as<double>();
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
	M_throw() << "Failed a lexical cast in CLDblWall";
      }
  }

  void 
  CLDblWall::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "DoubleWall" 
	<< magnet::xml::attr("Name") << localName
	<< magnet::xml::attr("Elasticity") << e
	<< range
	<< magnet::xml::tag("Norm")
	<< vNorm
	<< magnet::xml::endtag("Norm")
	<< magnet::xml::tag("Origin")
	<< vPosition / Sim->dynamics.units().unitLength()
	<< magnet::xml::endtag("Origin");
  }
}
