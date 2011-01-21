/*  DYNAMO:- Event driven molecular dynamics simulator 
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

#include "intEvent.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../../simulation/particle.hpp"
#include "../../base/is_simdata.hpp"
#include "../units/units.hpp"
#include <cmath>

xml::XmlStream& operator<<(xml::XmlStream &XML, 
			    const IntEvent &coll)
{
  XML << xml::tag("Collision")
      << xml::attr("p1ID") << coll.getParticle1ID()
      << coll.particle1
      << xml::attr("p2ID") << coll.getParticle1ID()
      << coll.particle2
      << xml::attr("dt")   << coll.dt
      << xml::endtag("Collision");
  
  return XML;
}

std::string 
IntEvent::stringData(const DYNAMO::SimData* Sim) const
{
  std::ostringstream tmpstring;
  tmpstring << "dt :" << dt / Sim->dynamics.units().unitTime()
	    << "\nType :" << CType
	    << "\nP1 :" << particle1;

  if (hasParticle2())
    tmpstring << "\nP2 :" << particle2;

  return tmpstring.str();
}
