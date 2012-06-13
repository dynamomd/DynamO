/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include <dynamo/globals/include.hpp>
#include <dynamo/globals/global.hpp>
#include <dynamo/particle.hpp>
#include <dynamo/globals/globEvent.hpp>
#include <dynamo/ranges/1RAll.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
  Global::Global(dynamo::Simulation* tmp, std::string name, Range* nR):
    SimBase(tmp, name),
    range(nR ? nR : new RAll(tmp))
  {}

  bool 
  Global::isInteraction(const Particle &p1) const
  {
    return range->isInRange(p1);
  }

  magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream& XML, const Global& g)
  {
    g.outputXML(XML);
    return XML;
  }

  shared_ptr<Global>
  Global::getClass(const magnet::xml::Node& XML, dynamo::Simulation* Sim)
  {
    if (!strcmp(XML.getAttribute("Type"),"Cells2")
	|| !strcmp(XML.getAttribute("Type"),"Cells")
	|| !strcmp(XML.getAttribute("Type"),"CellsMorton"))
      return shared_ptr<Global>(new GCells(XML, Sim));
    else if (!strcmp(XML.getAttribute("Type"),"ShearingCells"))
      return shared_ptr<Global>(new GCellsShearing(XML, Sim));
    else if (!strcmp(XML.getAttribute("Type"),"SOCells"))
      return shared_ptr<Global>(new GSOCells(XML, Sim));
    else if (!strcmp(XML.getAttribute("Type"),"Waker"))
      return shared_ptr<Global>(new GWaker(XML, Sim));
    else 
      M_throw() << XML.getAttribute("Type")
		<< ", Unknown type (" << XML.getAttribute("Type")
		<< ") of Global Interaction encountered";
  }
}

