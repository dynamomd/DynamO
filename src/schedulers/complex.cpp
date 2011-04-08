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

#include "complex.hpp"
#include "../dynamics/interactions/intEvent.hpp"
#include "../simulation/particle.hpp"
#include "../dynamics/dynamics.hpp"
#include "../dynamics/liouvillean/liouvillean.hpp"
#include "../base/is_simdata.hpp"
#include "../base/is_base.hpp"
#include "../dynamics/systems/system.hpp"
#include "../dynamics/globals/global.hpp"
#include "../dynamics/globals/globEvent.hpp"
#include "../dynamics/globals/neighbourList.hpp"
#include "../dynamics/locals/local.hpp"
#include "../dynamics/locals/localEvent.hpp"
#include "complexentries/include.hpp"

#include <magnet/xmlreader.hpp>
#include <boost/bind.hpp>
#include <boost/progress.hpp>
#include <cmath> //for huge val


void 
CSComplex::operator<<(const magnet::xml::Node& XML)
{
  sorter.set_ptr(CSSorter::getClass(XML.getNode("Sorter"), Sim));

  for (magnet::xml::Node node = XML.getNode("Entries").getNode("Entry"); 
       node.valid(); ++node)
    entries.push_back(magnet::ClonePtr<CSCEntry>(CSCEntry::getClass(node, Sim)));
}

void
CSComplex::initialise()
{
  I_cout() << "Reinitialising on collision " << Sim->eventCount;
  std::cout.flush();

  BOOST_FOREACH(magnet::ClonePtr<CSCEntry>& ent, entries)
    ent->initialise();

  sorter->clear();

  //The plus one is because system events are stored in the last heap;
  sorter->resize(Sim->N+1);
  eventCount.clear();
  eventCount.resize(Sim->N+1, 0);

  //Now initialise the interactions
  {
    boost::progress_display prog(Sim->N);
 
    BOOST_FOREACH(const Particle& part, Sim->particleList)
      {
	addEventsInit(part);
	++prog;
      }
  }
  
  sorter->init();

  rebuildSystemEvents();
}

void 
CSComplex::rebuildList()
{
#ifdef DYNAMO_DEBUG
  initialise();
#else
  BOOST_FOREACH(magnet::ClonePtr<CSCEntry>& ent, entries)
    ent->initialise();

  sorter->clear();

  //The plus one is because system events are stored in the last heap;
  sorter->resize(Sim->N+1);
  eventCount.clear();
  eventCount.resize(Sim->N+1, 0);

  BOOST_FOREACH(const Particle& part, Sim->particleList)
    addEventsInit(part);
  
  sorter->rebuild();

  rebuildSystemEvents();
#endif
}


void 
CSComplex::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Type") << "Complex"
      << xml::tag("Sorter")
      << sorter
      << xml::endtag("Sorter")
      << xml::tag("Entries");
  
  BOOST_FOREACH(const magnet::ClonePtr<CSCEntry>& ent,  entries)
    XML << xml::tag("Entry")
	<< ent
	<< xml::endtag("Entry");

  XML << xml::endtag("Entries");
}

CSComplex::CSComplex(const magnet::xml::Node& XML, DYNAMO::SimData* const Sim):
  CScheduler(Sim,"ComplexScheduler", NULL)
{ 
  I_cout() << "Complex Scheduler Algorithmn Loaded";
  operator<<(XML);
}

CSComplex::CSComplex(DYNAMO::SimData* const Sim, CSSorter* ns):
  CScheduler(Sim,"ComplexScheduler", ns)
{ I_cout() << "Complex Scheduler Algorithmn Loaded"; }

void 
CSComplex::addEvents(const Particle& part)
{
  Sim->dynamics.getLiouvillean().updateParticle(part);
  
  //Add the global events
  BOOST_FOREACH(const magnet::ClonePtr<Global>& glob, Sim->dynamics.getGlobals())
    if (glob->isInteraction(part))
      sorter->push(glob->getEvent(part), part.getID());
  
  BOOST_FOREACH(const magnet::ClonePtr<CSCEntry>& ent, entries)
    if (ent->isApplicable(part))
      {
	ent->getParticleLocalNeighbourhood
	  (part, magnet::function::MakeDelegate(this, &CScheduler::addLocalEvent));

	ent->getParticleNeighbourhood
	  (part, magnet::function::MakeDelegate(this, &CScheduler::addInteractionEvent));
      }
}

void 
CSComplex::addEventsInit(const Particle& part)
{  
  Sim->dynamics.getLiouvillean().updateParticle(part);

  //Add the global events
  BOOST_FOREACH(const magnet::ClonePtr<Global>& glob, Sim->dynamics.getGlobals())
    if (glob->isInteraction(part))
      sorter->push(glob->getEvent(part), part.getID());
  
  BOOST_FOREACH(const magnet::ClonePtr<CSCEntry>& ent, entries)
    if (ent->isApplicable(part))
      {
	ent->getParticleLocalNeighbourhood
	  (part, magnet::function::MakeDelegate(this, &CScheduler::addLocalEvent));

	ent->getParticleNeighbourhood
	  (part, magnet::function::MakeDelegate(this, &CScheduler::addInteractionEventInit));
      }
}
