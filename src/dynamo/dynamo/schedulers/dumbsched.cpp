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

#include <dynamo/schedulers/dumbsched.hpp>
#include <dynamo/interactions/intEvent.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/locals/local.hpp>
#include <dynamo/ranges/1RAll.hpp>
#include <magnet/xmlreader.hpp>
#include <cmath> //for huge val

namespace dynamo {
  SDumb::SDumb(const magnet::xml::Node& XML, dynamo::Simulation* const Sim):
    Scheduler(Sim,"DumbScheduler", NULL)
  { 
    dout << "Dumb Scheduler Algorithmn" << std::endl;
    operator<<(XML);
  }

  SDumb::SDumb(dynamo::Simulation* const Sim, FEL* ns):
    Scheduler(Sim,"DumbScheduler", ns)
  { dout << "Dumb Scheduler Algorithmn" << std::endl; }

  void 
  SDumb::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "Dumb"
	<< magnet::xml::tag("Sorter")
	<< *sorter
	<< magnet::xml::endtag("Sorter");
  }

  std::auto_ptr<Range>
  SDumb::getParticleNeighbours(const Particle&) const
  {
    return std::auto_ptr<Range>(new RAll(Sim));
  }

  void 
  SDumb::getParticleNeighbourhood(const Particle& part,
				  const nbHoodFunc& func) const
  {
    BOOST_FOREACH(const Particle& op, Sim->particles)
      if (op != part) func(part, op.getID());
  }
    
  void 
  SDumb::getLocalNeighbourhood(const Particle& part, 
				       const nbHoodFunc& func) const
  {
    BOOST_FOREACH(const shared_ptr<Local>& local, Sim->locals)
      func(part, local->getID());
  }

  void 
  SDumb::getParticleNeighbourhood(const Vector& vec, const nbHoodFunc2& func) const
  {
    BOOST_FOREACH(const Particle& op, Sim->particles)
      func(op.getID());
  }
}
