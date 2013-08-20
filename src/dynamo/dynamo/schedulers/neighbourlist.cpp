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

#include <dynamo/schedulers/neighbourlist.hpp>
#include <dynamo/interactions/intEvent.hpp>
#include <dynamo/particle.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/systems/system.hpp>
#include <dynamo/globals/global.hpp>
#include <dynamo/globals/globEvent.hpp>
#include <dynamo/globals/neighbourList.hpp>
#include <dynamo/locals/local.hpp>
#include <dynamo/locals/localEvent.hpp>
#include <dynamo/ranges/IDRangeRange.hpp>
#include <magnet/xmlreader.hpp>
#include <cmath>

namespace dynamo {
  void
  SNeighbourList::initialise()
  {
    try {
      NBListID = Sim->globals["SchedulerNBList"]->getID();
    }
    catch(std::exception& cxp)
      {
	M_throw() << "Failed while finding the neighbour list global.\n"
		  << "You must have a neighbour list enabled for this\n"
		  << "scheduler called SchedulerNBList.\n"
		  << cxp.what();
      }
  
    shared_ptr<GNeighbourList> nblist 
      = std::dynamic_pointer_cast<GNeighbourList>(Sim->globals[NBListID]);

    if (!nblist)
      M_throw() << "The Global named SchedulerNBList is not a neighbour list!";

    if (nblist->getMaxSupportedInteractionLength() 
	< Sim->getLongestInteraction())
      M_throw() << "Neighbourlist supports too small interaction distances! Supported distance is " 
		<< nblist->getMaxSupportedInteractionLength() / Sim->units.unitLength() 
		<< " but the longest interaction distance is " 
		<< Sim->getLongestInteraction() / Sim->units.unitLength();

    nblist->markAsUsedInScheduler();
    nblist->_sigNewNeighbour.connect<Scheduler, &Scheduler::addInteractionEvent>(this);
    Scheduler::initialise();
  }

  void 
  SNeighbourList::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "NeighbourList"
	<< magnet::xml::tag("Sorter")
	<< *sorter
	<< magnet::xml::endtag("Sorter");
  }

  SNeighbourList::SNeighbourList(const magnet::xml::Node& XML, 
				   dynamo::Simulation* const Sim):
    Scheduler(Sim,"NbListScheduler", NULL)
  { 
    dout << "Neighbour List Scheduler Algorithmn Loaded" << std::endl;
    operator<<(XML);
  }

  SNeighbourList::SNeighbourList(dynamo::Simulation* const Sim, FEL* ns):
    Scheduler(Sim,"NeighbourListScheduler", ns)
  { dout << "Neighbour List Scheduler Algorithmn Loaded" << std::endl; }

  std::unique_ptr<IDRange>
  SNeighbourList::getParticleNeighbours(const Particle& part) const
  {
#ifdef DYNAMO_DEBUG
    if (!std::dynamic_pointer_cast<GNeighbourList>(Sim->globals[NBListID]))
      M_throw() << "Not a GNeighbourList!";
#endif

    //Grab a reference to the neighbour list
    const GNeighbourList& nblist(*static_cast<const GNeighbourList*>
				 (Sim->globals[NBListID]
				  .get()));
  
    IDRangeList* range_ptr = new IDRangeList();
    nblist.getParticleNeighbours(part, range_ptr->getContainer());
    return std::unique_ptr<IDRange>(range_ptr);
  }

  std::unique_ptr<IDRange>
  SNeighbourList::getParticleNeighbours(const Vector& vec) const
  {
#ifdef DYNAMO_DEBUG
    if (!std::dynamic_pointer_cast<GNeighbourList>(Sim->globals[NBListID]))
      M_throw() << "Not a GNeighbourList!";
#endif

    //Grab a reference to the neighbour list
    const GNeighbourList& nblist(*static_cast<const GNeighbourList*>
				 (Sim->globals[NBListID]
				  .get()));
  
    IDRangeList* range_ptr = new IDRangeList();
    nblist.getParticleNeighbours(vec, range_ptr->getContainer());
    return std::unique_ptr<IDRange>(range_ptr);
  }
    
  std::unique_ptr<IDRange> 
  SNeighbourList::getParticleLocals(const Particle& part) const {
    return std::unique_ptr<IDRange>(new IDRangeRange(0, Sim->locals.size()));
  }
}
