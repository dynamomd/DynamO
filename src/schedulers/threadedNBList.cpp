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

#include "threadedNBList.hpp"
#include "../dynamics/interactions/intEvent.hpp"
#include "../simulation/particle.hpp"
#include "../dynamics/dynamics.hpp"
#include "../dynamics/liouvillean/liouvillean.hpp"
#include "../base/is_simdata.hpp"
#include "../base/is_base.hpp"
#include "../dynamics/systems/system.hpp"
#include <cmath> //for huge val
#include "../extcode/xmlParser.h"
#include "../dynamics/globals/global.hpp"
#include "../dynamics/globals/globEvent.hpp"
#include "../dynamics/globals/neighbourList.hpp"
#include "../dynamics/locals/local.hpp"
#include "../dynamics/locals/localEvent.hpp"
#include <boost/bind.hpp>
#include <boost/progress.hpp>

SThreadedNBList::SThreadedNBList(const XMLNode& XML, 
				 DYNAMO::SimData* const Sim):
  CSNeighbourList(XML, Sim)
{ 
  I_cout() << "Threaded Variant Loaded";
}

SThreadedNBList::SThreadedNBList(DYNAMO::SimData* const Sim, CSSorter* ns):
  CSNeighbourList(Sim, ns)
{ I_cout() << "Threaded Variant Loaded"; }


void 
SThreadedNBList::operator<<(const XMLNode& XML)
{
  CSNeighbourList::operator<<(XML);
}

void 
SThreadedNBList::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "ThreadedNeighbourList"
      << xmlw::tag("Sorter")
      << sorter
      << xmlw::endtag("Sorter");
}

void 
SThreadedNBList::addEvents(const CParticle& part)
{
  Sim->dynamics.getLiouvillean().updateParticle(part);
  
  //Add the global events
  BOOST_FOREACH(const smrtPlugPtr<CGlobal>& glob, Sim->dynamics.getGlobals())
    if (glob->isInteraction(part))
      sorter->push(glob->getEvent(part), part.getID());
  
#ifdef DYNAMO_DEBUG
  if (dynamic_cast<const CGNeighbourList*>
      (Sim->dynamics.getGlobals()[NBListID].get_ptr())
      == NULL)
    D_throw() << "Not a CGNeighbourList!";
#endif

  //Grab a reference to the neighbour list
  const CGNeighbourList& nblist(*static_cast<const CGNeighbourList*>
				(Sim->dynamics.getGlobals()[NBListID]
				 .get_ptr()));
  
  //Add the local cell events
  nblist.getParticleLocalNeighbourhood
    (part, fastdelegate::MakeDelegate(this, &CScheduler::addLocalEvent));

  //Add the interaction events
  nblist.getParticleNeighbourhood
    (part, fastdelegate::MakeDelegate(this, &SThreadedNBList::streamParticles));  

  nblist.getParticleNeighbourhood
    (part, fastdelegate::MakeDelegate(this, &SThreadedNBList::addEvents2));  
}

void 
SThreadedNBList::streamParticles(const CParticle& part, 
				 const size_t& id) const
{
  Sim->dynamics.getLiouvillean().updateParticle(Sim->vParticleList[id]);
}

void 
SThreadedNBList::addEvents2(const CParticle& part, 
			    const size_t& id) const
{
  const CIntEvent& eevent(Sim->dynamics.getEvent(part, Sim->vParticleList[id]));
  
  if (eevent.getType() != NONE)
    sorter->push(intPart(eevent, eventCount[id]), part.getID());
}


void 
SThreadedNBList::addEventsInit(const CParticle& part)
{  
  Sim->dynamics.getLiouvillean().updateParticle(part);

  //Add the global events
  BOOST_FOREACH(const smrtPlugPtr<CGlobal>& glob, Sim->dynamics.getGlobals())
    if (glob->isInteraction(part))
      sorter->push(glob->getEvent(part), part.getID());
  
#ifdef DYNAMO_DEBUG
  if (dynamic_cast<const CGNeighbourList*>
      (Sim->dynamics.getGlobals()[NBListID].get_ptr())
      == NULL)
    D_throw() << "Not a CGNeighbourList!";
#endif

  //Grab a reference to the neighbour list
  const CGNeighbourList& nblist(*static_cast<const CGNeighbourList*>
				(Sim->dynamics.getGlobals()[NBListID]
				 .get_ptr()));
  
  //Add the local cell events
  nblist.getParticleLocalNeighbourhood
    (part, fastdelegate::MakeDelegate(this, &CScheduler::addLocalEvent));

  //Add the interaction events
  nblist.getParticleNeighbourhood
    (part, fastdelegate::MakeDelegate(this, &CScheduler::addInteractionEventInit));  
}
