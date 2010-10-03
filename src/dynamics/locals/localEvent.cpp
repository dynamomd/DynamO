/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
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

#include "localEvent.hpp"
#include "local.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../../simulation/particle.hpp"
#include "../../base/is_simdata.hpp"
#include "../units/units.hpp"
#include "../interactions/intEvent.hpp"
#include <cmath>

LocalEvent::LocalEvent(const Particle& part1, const double &delt, 
			 EEventType nType, const Local& local):
  particle_(&part1), dt(delt), 
  CType(nType), localID(local.getID())
{}
  

xml::XmlStream& operator<<(xml::XmlStream &XML, 
			    const LocalEvent &coll)
{
  XML << xml::tag("Collision")
      << xml::attr("p1ID") << coll.getParticle().getID()
      << xml::attr("dt")   << coll.dt
      << xml::endtag("Collision");
  
  return XML;
}

std::string 
LocalEvent::stringData(const DYNAMO::SimData* Sim) const
{
  std::ostringstream tmpstring;
  tmpstring << "dt :" << dt / Sim->dynamics.units().unitTime()
	    << "\nType :" << CType
	    << "\nP1 :" << particle_->getID();
    return tmpstring.str();
}

bool 
LocalEvent::areInvolved(const IntEvent &coll) const 
{ return (coll == *particle_); }
