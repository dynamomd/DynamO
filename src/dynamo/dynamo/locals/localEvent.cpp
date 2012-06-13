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

#include <dynamo/locals/localEvent.hpp>
#include <dynamo/locals/local.hpp>
#include <dynamo/particle.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/interactions/intEvent.hpp>
#include <magnet/xmlwriter.hpp>
#include <cmath>

namespace dynamo {
  LocalEvent::LocalEvent(const Particle& part1, const double &delt, 
			 EEventType nType, const Local& local,
			 const size_t extraData):
    particle_(&part1), dt(delt), 
    CType(nType), localID(local.getID()),
    _extraData(extraData)
  {}
  

  magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream &XML, 
				     const LocalEvent &coll)
  {
    XML << magnet::xml::tag("Collision")
	<< magnet::xml::attr("p1ID") << coll.getParticle().getID()
	<< magnet::xml::attr("dt")   << coll.dt
	<< magnet::xml::endtag("Collision");
  
    return XML;
  }

  std::string 
  LocalEvent::stringData(const dynamo::Simulation* Sim) const
  {
    std::ostringstream tmpstring;
    tmpstring << "dt :" << dt / Sim->units.unitTime()
	      << "\nType :" << CType
	      << "\nP1 :" << particle_->getID();
    return tmpstring.str();
  }

  bool 
  LocalEvent::areInvolved(const IntEvent &coll) const 
  { return (coll == *particle_); }
}

