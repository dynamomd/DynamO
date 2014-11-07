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
#include <dynamo/locals/boundary/object.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
  namespace boundary {
    Object::Object(Simulation* const SD, const std::string aName):
      dynamo::SimBase(SD, aName)
    {}

    shared_ptr<Object>
    Object::getClass(const magnet::xml::Node& XML, dynamo::Simulation* Sim) {
      if (!XML.getAttribute("Type").getValue().compare("PlanarWall"))
	return shared_ptr<Object>(new PlanarWall(XML, Sim));
      else
	M_throw() << XML.getAttribute("Type").getValue()
		  << ", Unknown type of Object encountered" 
		  << XML.getPath();
    }

    PlanarWall::PlanarWall(const magnet::xml::Node& XML, dynamo::Simulation* Sim):
      Object(Sim, "PlanarWall")
    {}

    void
    PlanarWall::outputXML(magnet::xml::XmlStream& XML) const {
      XML << magnet::xml::attr("Type") << "PlanarWall";
    }

    bool
    PlanarWall::validateState(const Particle& part, bool textoutput) const {
      return true;
    }
    
    Event 
    PlanarWall::getEvent(const Particle& part) const {
      return Event();
    }
    
    ParticleEventData 
    PlanarWall::runEvent(Particle& part, const Event& event) const {
      return ParticleEventData();
    }
  }
}
