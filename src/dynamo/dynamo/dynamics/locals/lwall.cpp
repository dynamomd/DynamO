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

#include <dynamo/dynamics/locals/lwall.hpp>
#include <dynamo/dynamics/BC/BC.hpp>
#include <dynamo/dynamics/liouvillean/liouvillean.hpp>
#include <dynamo/dynamics/locals/localEvent.hpp>
#include <dynamo/dynamics/NparticleEventData.hpp>
#include <dynamo/dynamics/units/units.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <magnet/overlap/cube_plane.hpp>

namespace dynamo {
  LWall::LWall(const magnet::xml::Node& XML, dynamo::SimData* tmp):
    Local(tmp, "LocalWall")
  {
    operator<<(XML);
  }

  LocalEvent 
  LWall::getEvent(const Particle& part) const
  {
#ifdef ISSS_DEBUG
    if (!Sim->dynamics.getLiouvillean().isUpToDate(part))
      M_throw() << "Particle is not up to date";
#endif

    double colldist = 0.5 * _diameter->getProperty(part.getID());

    return LocalEvent(part, Sim->dynamics.getLiouvillean().getWallCollision
		      (part, vPosition + vNorm * colldist, vNorm), WALL, *this);
  }

  void
  LWall::runEvent(const Particle& part, const LocalEvent& iEvent) const
  {
    ++Sim->eventCount;

    //Run the collision and catch the data

    NEventData EDat(Sim->dynamics.getLiouvillean().runWallCollision
		    (part, vNorm, _e->getProperty(part.getID())));

    Sim->signalParticleUpdate(EDat);

    //Now we're past the event update the scheduler and plugins
    Sim->ptrScheduler->fullUpdate(part);
  
    BOOST_FOREACH(shared_ptr<OutputPlugin> & Ptr, Sim->outputPlugins)
      Ptr->eventUpdate(iEvent, EDat);
  }

  bool 
  LWall::isInCell(const Vector & Origin, const Vector& CellDim) const
  {
    double max_int_dist = 0.5 * _diameter->getMaxValue();

    Vector data(max_int_dist, max_int_dist, max_int_dist);

    return magnet::overlap::cube_plane(Origin - 0.5 * data, CellDim + data, 
				       vPosition, vNorm);
  }

  void 
  LWall::initialise(size_t nID)
  {
    ID = nID;
  }

  void 
  LWall::operator<<(const magnet::xml::Node& XML)
  {
    range = shared_ptr<Range>(Range::getClass(XML,Sim));
  
    try {
      _diameter = Sim->_properties.getProperty(XML.getAttribute("Diameter"),
					       Property::Units::Length());
      _e = Sim->_properties.getProperty(XML.getAttribute("Elasticity"),
					Property::Units::Dimensionless());
    
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
	M_throw() << "Failed a lexical cast in LWall";
      }
  }

  void 
  LWall::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "Wall" 
	<< magnet::xml::attr("Name") << localName
	<< magnet::xml::attr("Elasticity") << _e->getName()
	<< magnet::xml::attr("Diameter") << _diameter->getName()
	<< range
	<< magnet::xml::tag("Norm")
	<< vNorm
	<< magnet::xml::endtag("Norm")
	<< magnet::xml::tag("Origin")
	<< vPosition / Sim->dynamics.units().unitLength()
	<< magnet::xml::endtag("Origin");
  }

  void 
  LWall::checkOverlaps(const Particle& p1) const
  {
    Vector pos(p1.getPosition() - vPosition);
    Sim->dynamics.BCs().applyBC(pos);

    double r = (pos | vNorm);
  
    if (r < 0)
      dout << "Possible overlap of " << r / Sim->dynamics.units().unitLength() << " for particle " << p1.getID()
	   << "\nWall Pos is [" 
	   << vPosition[0] << "," << vPosition[1] << "," << vPosition[2] 
	   << "] and Normal is [" 
	   << vNorm[0] << "," << vNorm[1] << "," << vNorm[2] << "]"
	   << std::endl;
  }
}

