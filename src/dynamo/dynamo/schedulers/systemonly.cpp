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

#include "systemonly.hpp"
#include "../dynamics/interactions/intEvent.hpp"
#include "../simulation/particle.hpp"
#include "../dynamics/dynamics.hpp"
#include "../dynamics/liouvillean/liouvillean.hpp"
#include "../dynamics/BC/BC.hpp"
#include "../dynamics/BC/LEBC.hpp"
#include "../base/is_simdata.hpp"
#include "../dynamics/globals/globEvent.hpp"
#include "../dynamics/systems/system.hpp"
#include "../dynamics/globals/global.hpp"
#include "../dynamics/globals/globEvent.hpp"
#include "../dynamics/globals/neighbourList.hpp"
#include "../dynamics/locals/local.hpp"
#include "../dynamics/locals/localEvent.hpp"
#include <magnet/xmlreader.hpp>
#include <cmath> //for huge val

CSSystemOnly::CSSystemOnly(const magnet::xml::Node& XML, 
			   dynamo::SimData* const Sim):
  CScheduler(Sim,"SystemOnlyScheduler", NULL)
{ 
  dout << "System Events Only Scheduler Algorithmn" << std::endl;
  operator<<(XML);
}

CSSystemOnly::CSSystemOnly(dynamo::SimData* const Sim, CSSorter* ns):
  CScheduler(Sim,"SystemOnlyScheduler", ns)
{ dout << "System Events Only Scheduler Algorithmn" << std::endl; }

void 
CSSystemOnly::operator<<(const magnet::xml::Node& XML)
{
  sorter.set_ptr(CSSorter::getClass(XML.getNode("Sorter"), Sim));
}

void
CSSystemOnly::initialise()
{
  dout << "Reinitialising on collision " << Sim->eventCount << std::endl;

  if (Sim->dynamics.getSystemEvents().empty())
    M_throw() << "A SystemOnlyScheduler used when there are no system events?";
  
  sorter->clear();
  sorter->resize(Sim->N+1);
  eventCount.clear();
  eventCount.resize(Sim->N+1, 0);  
  sorter->init();
  rebuildSystemEvents();
}

void
CSSystemOnly::rebuildList()
{
#ifdef DYNAMO_DEBUG
  initialise();
#else
  if (Sim->dynamics.getSystemEvents().empty())
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
CSSystemOnly::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Type") << "SystemOnly"
      << xml::tag("Sorter")
      << sorter
      << xml::endtag("Sorter");
}

void 
CSSystemOnly::addEvents(const Particle& part)
{}
