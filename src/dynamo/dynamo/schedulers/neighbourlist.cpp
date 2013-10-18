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
#include <dynamo/dynamics/compression.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/systems/system.hpp>
#include <dynamo/globals/cells.hpp>
#include <dynamo/globals/cellsShearing.hpp>
#include <dynamo/systems/nblistCompressionFix.hpp>
#include <dynamo/locals/local.hpp>
#include <dynamo/locals/localEvent.hpp>
#include <dynamo/ranges/IDRangeRange.hpp>
#include <dynamo/BC/include.hpp>
#include <magnet/xmlreader.hpp>
#include <cmath>

namespace dynamo {
  void
  SNeighbourList::initialise()
  {
    //First, try to detect a neighbour list
    try {
      NBListID = Sim->globals["SchedulerNBList"]->getID();
    }
    catch(std::exception& cxp)
      {
	//There is no global cellular list available. Add the appropriate neighbourlist.
	shared_ptr<GCells> scheduler;
	if (std::dynamic_pointer_cast<BCLeesEdwards>(Sim->BCs))
	  scheduler = shared_ptr<GCells>(new GCellsShearing(Sim, "SchedulerNBList"));
	else
	  scheduler = shared_ptr<GCells>(new GCells(Sim,"SchedulerNBList"));
	Sim->globals.push_back(scheduler);
	scheduler->setConfigOutput(false);
	scheduler->initialise(Sim->globals.size() - 1);

	NBListID = Sim->globals.back()->getID();

	//Check if this is a compressing system, if so, add the appropriate hack!
	shared_ptr<DynCompression> compressiondynamics = std::dynamic_pointer_cast<DynCompression>(Sim->dynamics);
	if (compressiondynamics)
	  {
	    //Rebulid the collision scheduler without the overlapping
	    //cells, otherwise cells are always rebuilt as they overlap
	    //such that the maximum supported interaction distance is
	    //equal to the current maximum interaction distance.
	    scheduler->setCellOverlap(false);
	    Sim->systems.push_back(shared_ptr<System>(new SysNBListCompressionFix(Sim, compressiondynamics->getGrowthRate(), NBListID)));
	    Sim->systems.back()->initialise(Sim->systems.size() - 1);
	  }

      }
    
    shared_ptr<GNeighbourList> nblist = std::dynamic_pointer_cast<GNeighbourList>(Sim->globals[NBListID]);

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
    dout << "Neighbour List Scheduler Algorithm Loaded" << std::endl;
    operator<<(XML);
  }

  SNeighbourList::SNeighbourList(dynamo::Simulation* const Sim, FEL* ns):
    Scheduler(Sim,"NeighbourListScheduler", ns)
  { dout << "Neighbour List Scheduler Algorithm Loaded" << std::endl; }

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
    return std::unique_ptr<IDRange>(new IDRangeRange(0, Sim->locals.size() - 1));
  }
}
