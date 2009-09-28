/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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
#include "../base/is_base.hpp"
#include "../dynamics/globals/globEvent.hpp"
#include "../dynamics/systems/system.hpp"
#include <cmath> //for huge val
#include "../extcode/xmlParser.h"
#include "../dynamics/globals/global.hpp"
#include "../dynamics/globals/globEvent.hpp"
#include "../dynamics/globals/neighbourList.hpp"
#include "../dynamics/locals/local.hpp"
#include "../dynamics/locals/localEvent.hpp"

CSSystemOnly::CSSystemOnly(const XMLNode& XML, 
				 DYNAMO::SimData* const Sim):
  CScheduler(Sim,"SystemOnlyScheduler", NULL)
{ 
  I_cout() << "System Events Only Scheduler Algorithmn";
  operator<<(XML);
}

CSSystemOnly::CSSystemOnly(DYNAMO::SimData* const Sim, CSSorter* ns):
  CScheduler(Sim,"SystemOnlyScheduler", ns)
{ I_cout() << "System Events Only Scheduler Algorithmn"; }

void 
CSSystemOnly::operator<<(const XMLNode& XML)
{
  sorter.set_ptr(CSSorter::getClass(XML.getChildNode("Sorter"), Sim));
}

void
CSSystemOnly::initialise()
{
  I_cout() << "Reinitialising on collision " << Sim->lNColl;

  if (Sim->Dynamics.getSystemEvents().empty())
    D_throw() << "A SystemOnlyScheduler used when there are no system events?";
  
  sorter->clear();
  sorter->resize(Sim->lN+1);
  eventCount.clear();
  eventCount.resize(Sim->lN+1, 0);  
  sorter->init();
  rebuildSystemEvents();

}

void 
CSSystemOnly::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "SystemOnly"
      << xmlw::tag("Sorter")
      << sorter
      << xmlw::endtag("Sorter");
}

void 
CSSystemOnly::addEvents(const CParticle& part)
{}
