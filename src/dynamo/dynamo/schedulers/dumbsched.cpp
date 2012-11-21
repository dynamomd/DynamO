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
#include <dynamo/ranges/IDRangeAll.hpp>
#include <dynamo/ranges/IDRangeRange.hpp>
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

  std::auto_ptr<IDRange>
  SDumb::getParticleNeighbours(const Particle&) const
  {
    return std::auto_ptr<IDRange>(new IDRangeAll(Sim));
  }

  std::auto_ptr<IDRange>
  SDumb::getParticleNeighbours(const Vector&) const
  {
    return std::auto_ptr<IDRange>(new IDRangeAll(Sim));
  }

  std::auto_ptr<IDRange>
  SDumb::getParticleLocals(const Particle&) const
  {
    return std::auto_ptr<IDRange>(new IDRangeRange(0, Sim->locals.size()));
  }
}
