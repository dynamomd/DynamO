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

#include <dynamo/locals/lwall.hpp>
#include <dynamo/locals/oscillatingplate.hpp>
#include <dynamo/locals/lcylinder.hpp>
#include <dynamo/locals/lroughwall.hpp>
#include <dynamo/locals/trianglemesh.hpp>
#include <dynamo/locals/localEvent.hpp>
#include <dynamo/ranges/1RAll.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
  Local::Local(dynamo::Simulation* tmp, const char *name):
    SimBase(tmp,name),
    range(new RAll(tmp))
  {}

  Local::Local(Range* nR, dynamo::Simulation* tmp, const char *name):
    SimBase(tmp, name),
    range(nR)
  {}

  bool 
  Local::isInteraction(const Particle &p1) const
  {
    return range->isInRange(p1);
  }

  magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream& XML, 
				     const Local& g)
  {
    g.outputXML(XML);
    return XML;
  }

  shared_ptr<Local> 
  Local::getClass(const magnet::xml::Node& XML, dynamo::Simulation* Sim)
  {
    if (!strcmp(XML.getAttribute("Type"),"Wall"))
      return shared_ptr<Local>(new LWall(XML, Sim));
    else if (!strcmp(XML.getAttribute("Type"),"RoughWall"))
      return shared_ptr<Local>(new LRoughWall(XML, Sim));
    else if (!strcmp(XML.getAttribute("Type"),"TriangleMesh"))
      return shared_ptr<Local>(new LTriangleMesh(XML, Sim));
    else if (!strcmp(XML.getAttribute("Type"),"OscillatingPlate"))
      return shared_ptr<Local>(new LOscillatingPlate(XML, Sim));
    else if (!strcmp(XML.getAttribute("Type"),"CylinderWall"))
      return shared_ptr<Local>(new LCylinder(XML, Sim));
    else 
      M_throw() << XML.getAttribute("Type")
		<< ", Unknown type of Local Interaction encountered";
  }
}
