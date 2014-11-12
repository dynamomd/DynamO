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
    Object::Object(Simulation* const SD, const std::string aName, const BoundaryOscillationData& data):
      dynamo::SimBase(SD, aName),
      _oscillationData(data)
    {}

    shared_ptr<Object>
    Object::getClass(const magnet::xml::Node& XML, dynamo::Simulation* Sim, const BoundaryOscillationData& data) {
      if (!XML.getAttribute("Type").getValue().compare("PlanarWall"))
	return shared_ptr<Object>(new PlanarWall(XML, Sim, data));
      else
	M_throw() << XML.getAttribute("Type").getValue()
		  << ", Unknown type of Object encountered" 
		  << XML.getPath();
    }

    PlanarWall::PlanarWall(const magnet::xml::Node& XML, dynamo::Simulation* Sim, const BoundaryOscillationData& data):
      Object(Sim, "PlanarWall", data)
    {
      if (XML.hasNode("Position"))
	_position << XML.getNode("Position");
      _position *= Sim->units.unitLength();
      _normal << XML.getNode("Normal");      
    }

    void
    PlanarWall::outputXML(magnet::xml::XmlStream& XML) const {
      XML << magnet::xml::attr("Type") << "PlanarWall";
      if (_position.nrm() != 0)
	XML << magnet::xml::tag("Position")
	    << _position / Sim->units.unitLength()
	    << magnet::xml::endtag("Position");

      XML << magnet::xml::tag("Normal")
	  << _normal
	  << magnet::xml::endtag("Normal");
    }

    bool
    PlanarWall::validateState(const Particle& part, bool textoutput) const {
      return true;
    }
    
    Event 
    PlanarWall::getEvent(const Particle& part, const double diameter) const 
    {
      return Event();
    }
      
    Vector 
    PlanarWall::getContactNormal(const Particle&, const Event&) const
    {
      return _normal;
    }

    std::pair<std::vector<float>, std::vector<GLuint> >
    PlanarWall::getTessalatedSurfaces() const {
      const std::vector<float> vertices = {0,0,0,1, 1,0,0,0, 0,1,0,0, -1,0,0,0, 0,-1,0,0};
      const std::vector<GLuint> indices = {0,1,2, 0,2,3, 0,3,4, 0,4,1};
      return std::make_pair(vertices, indices);
    }
  }
}
