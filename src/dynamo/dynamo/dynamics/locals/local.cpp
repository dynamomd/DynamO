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

#include <dynamo/dynamics/locals/include.hpp>
#include <dynamo/dynamics/locals/local.hpp>
#include <dynamo/simulation/particle.hpp>
#include <dynamo/dynamics/locals/localEvent.hpp>
#include <dynamo/dynamics/ranges/1RAll.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
  Local::Local(dynamo::SimData* tmp, const char *name):
    SimBase(tmp,name),
    range(new RAll(tmp))
  {}

  Local::Local(Range* nR, dynamo::SimData* tmp, const char *name):
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
  Local::getClass(const magnet::xml::Node& XML, dynamo::SimData* Sim)
  {
    if (!strcmp(XML.getAttribute("Type"),"Wall"))
      return shared_ptr<Local>(new LWall(XML, Sim));
    else if (!strcmp(XML.getAttribute("Type"),"RoughWall"))
      return shared_ptr<Local>(new LRoughWall(XML, Sim));
    else if (!strcmp(XML.getAttribute("Type"),"TriangleMesh"))
      return shared_ptr<Local>(new LTriangleMesh(XML, Sim));
    else if (!strcmp(XML.getAttribute("Type"),"AndersenWall"))
      return shared_ptr<Local>(new LAndersenWall(XML, Sim));
    else if (!strcmp(XML.getAttribute("Type"),"DoubleWall"))
      return shared_ptr<Local>(new LDblWall(XML, Sim));
    else if (!strcmp(XML.getAttribute("Type"),"OscillatingPlate"))
      return shared_ptr<Local>(new LOscillatingPlate(XML, Sim));
    else if (!strcmp(XML.getAttribute("Type"),"CylinderWall"))
      return shared_ptr<Local>(new LCylinder(XML, Sim));
    else if (!strcmp(XML.getAttribute("Type"),"SphereWall"))
      return shared_ptr<Local>(new LSphere(XML, Sim));
    else 
      M_throw() << XML.getAttribute("Type")
		<< ", Unknown type of Local Interaction encountered";
  }
}
