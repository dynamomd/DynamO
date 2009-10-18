/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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
#include "../../base/is_exception.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../../simulation/particle.hpp"
#include "../../base/is_simdata.hpp"
#include "../units/units.hpp"
#include <cmath>

const char * 
CIntEvent::getCollEnumName(EEventType a)
{
  switch (a)
    {
    case NONE:        return "NONE";
    case CELL:        return "CELL";
    case GLOBAL:      return "GLOBAL";
    case INTERACTION: return "INTERACTION";
    case SYSTEM:      return "SYSTEM";
    case LOCAL:       return "LOCAL";
    case CORE:        return "CORE";
    case WELL_IN:     return "WELL_IN";
    case WELL_OUT:    return "WELL_OUT";
    case WELL_KEUP:   return "WELL_KEUP";
    case WELL_KEDOWN: return "WELL_KEDOWN";
    case BOUNCE:      return "BOUNCE";
    case WALL:        return "WALL";
    case GAUSSIAN:    return "GAUSSIAN";
    case DSMC:        return "DSMC";
    case HALT:        return "HALT";
    case STREAM:      return "STREAM";
    case NON_EVENT:   return "NON_EVENT";
    case VIRTUAL:     return "VIRTUAL";
    default:
      D_throw() << "Not a defined collision enumeration " << 
	a;
    }
}

xmlw::XmlStream& operator<<(xmlw::XmlStream &XML, 
			    const CIntEvent &coll)
{
  XML << xmlw::tag("Collision")
      << xmlw::attr("p1ID") << coll.getParticle1ID()
      << coll.particle1
      << xmlw::attr("p2ID") << coll.getParticle1ID()
      << coll.particle2
      << xmlw::attr("dt")   << coll.dt
      << xmlw::endtag("Collision");
  
  return XML;
}

std::string 
CIntEvent::stringData(const DYNAMO::SimData* Sim) const
{
  std::ostringstream tmpstring;
  tmpstring << "dt :" << dt / Sim->dynamics.units().unitTime()
	    << "\nType :" << getCollEnumName(CType)
	    << "\nP1 :" << particle1;

  if (hasParticle2())
    tmpstring << "\nP2 :" << particle2;

  return tmpstring.str();
}
