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

#include <dynamo/schedulers/systemonly.hpp>
#include <dynamo/interactions/intEvent.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/ranges/1RNone.hpp>
#include <magnet/xmlreader.hpp>
#include <cmath> //for huge val

namespace dynamo {
  SSystemOnly::SSystemOnly(const magnet::xml::Node& XML, 
			     dynamo::Simulation* const Sim):
    Scheduler(Sim,"SystemOnlyScheduler", NULL)
  { 
    dout << "System Events Only Scheduler Algorithmn" << std::endl;
    operator<<(XML);
  }

  SSystemOnly::SSystemOnly(dynamo::Simulation* const Sim, FEL* ns):
    Scheduler(Sim,"SystemOnlyScheduler", ns)
  { dout << "System Events Only Scheduler Algorithmn" << std::endl; }

  void
  SSystemOnly::initialise()
  {
    dout << "Reinitialising on collision " << Sim->eventCount << std::endl;

    if (Sim->systems.empty())
      M_throw() << "A SystemOnlyScheduler used when there are no system events?";
  
    sorter->clear();
    sorter->resize(Sim->N+1);
    eventCount.clear();
    eventCount.resize(Sim->N+1, 0);  
    sorter->init();
    rebuildSystemEvents();
  }

  void
  SSystemOnly::rebuildList()
  {
#ifdef DYNAMO_DEBUG
    initialise();
#else
    if (Sim->systems.empty())
      M_throw() << "A SystemOnlyScheduler used when there are no system events?";
  
    sorter->clear();
    sorter->resize(Sim->N+1);
    eventCount.clear();
    eventCount.resize(Sim->N+1, 0);  
    sorter->rebuild();
    rebuildSystemEvents();
#endif
  }

  void 
  SSystemOnly::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "SystemOnly"
	<< magnet::xml::tag("Sorter")
	<< *sorter
	<< magnet::xml::endtag("Sorter");
  }

  std::auto_ptr<Range>
  SSystemOnly::getParticleNeighbours(const Particle&) const
  {
    return std::auto_ptr<Range>(new RNone());
  }
}
