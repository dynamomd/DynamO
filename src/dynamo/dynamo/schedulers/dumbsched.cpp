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

#include <dynamo/schedulers/dumbsched.hpp>
#include <dynamo/dynamics/interactions/intEvent.hpp>
#include <dynamo/simulation/particle.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/dynamics/liouvillean/liouvillean.hpp>
#include <dynamo/dynamics/BC/BC.hpp>
#include <dynamo/dynamics/BC/LEBC.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <dynamo/dynamics/globals/globEvent.hpp>
#include <dynamo/dynamics/systems/system.hpp>
#include <dynamo/dynamics/globals/global.hpp>
#include <dynamo/dynamics/globals/globEvent.hpp>
#include <dynamo/dynamics/globals/neighbourList.hpp>
#include <dynamo/dynamics/locals/local.hpp>
#include <dynamo/dynamics/locals/localEvent.hpp>
#include <magnet/xmlreader.hpp>
#include <cmath> //for huge val

namespace dynamo {
  SDumb::SDumb(const magnet::xml::Node& XML, dynamo::SimData* const Sim):
    Scheduler(Sim,"DumbScheduler", NULL)
  { 
    dout << "Dumb Scheduler Algorithmn" << std::endl;
    operator<<(XML);
  }

  SDumb::SDumb(dynamo::SimData* const Sim, CSSorter* ns):
    Scheduler(Sim,"DumbScheduler", ns)
  { dout << "Dumb Scheduler Algorithmn" << std::endl; }

  void 
  SDumb::operator<<(const magnet::xml::Node& XML)
  {
    sorter.set_ptr(CSSorter::getClass(XML.getNode("Sorter"), Sim));
  }

  void
  SDumb::initialise()
  {
    dout << "Reinitialising on collision " << Sim->eventCount << std::endl;
  
    sorter->clear();
    sorter->resize(Sim->N+1);
    eventCount.clear();
    eventCount.resize(Sim->N+1, 0);
  
    //Now initialise the interactions
    BOOST_FOREACH(const Particle& part, Sim->particleList)
      addEvents(part);
  
    sorter->init();
    rebuildSystemEvents();
  }

  void 
  SDumb::rebuildList()
  {
#ifdef DYNAMO_DEBUG
    initialise();
#else
    sorter->clear();
    sorter->resize(Sim->N+1);
    eventCount.clear();
    eventCount.resize(Sim->N+1, 0);
  
    //Now initialise the interactions
    BOOST_FOREACH(const Particle& part, Sim->particleList)
      addEvents(part);
  
    sorter->rebuild();
    rebuildSystemEvents();
#endif
  }


  void 
  SDumb::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "Dumb"
	<< magnet::xml::tag("Sorter")
	<< sorter
	<< magnet::xml::endtag("Sorter");
  }

  void 
  SDumb::addEvents(const Particle& part)
  {  
    Sim->dynamics.getLiouvillean().updateParticle(part);

    //Add the global events
    BOOST_FOREACH(const magnet::ClonePtr<Global>& glob, Sim->dynamics.getGlobals())
      if (glob->isInteraction(part))
	sorter->push(glob->getEvent(part), part.getID());
  
    //Add the local cell events
    BOOST_FOREACH(const magnet::ClonePtr<Local>& local, Sim->dynamics.getLocals())
      addLocalEvent(part, local->getID());

    //Add the interaction events
    BOOST_FOREACH(const Particle& part2, Sim->particleList)
      if (part2 != part)
	addInteractionEvent(part, part2.getID());
  }
}
