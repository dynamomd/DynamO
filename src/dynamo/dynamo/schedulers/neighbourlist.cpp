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

#include <dynamo/schedulers/neighbourlist.hpp>
#include <dynamo/dynamics/interactions/intEvent.hpp>
#include <dynamo/simulation/particle.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/dynamics/liouvillean/liouvillean.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <dynamo/dynamics/systems/system.hpp>
#include <dynamo/dynamics/globals/global.hpp>
#include <dynamo/dynamics/globals/globEvent.hpp>
#include <dynamo/dynamics/globals/neighbourList.hpp>
#include <dynamo/dynamics/locals/local.hpp>
#include <dynamo/dynamics/locals/localEvent.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/bind.hpp>
#include <boost/progress.hpp>
#include <cmath>

namespace dynamo {
  void
  SNeighbourList::initialise()
  {
    try {
      NBListID = Sim->dynamics.getGlobal("SchedulerNBList")->getID();
    }
    catch(std::exception& cxp)
      {
	M_throw() << "Failed while finding the neighbour list global.\n"
		  << "You must have a neighbour list enabled for this\n"
		  << "scheduler called SchedulerNBList.\n"
		  << cxp.what();
      }
  
    std::tr1::shared_ptr<GNeighbourList> nblist 
      = std::tr1::dynamic_pointer_cast<GNeighbourList>(Sim->dynamics.getGlobals()[NBListID]);

    if (!nblist)
      M_throw() << "The Global named SchedulerNBList is not a neighbour list!";

    nblist->markAsUsedInScheduler();
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
				   dynamo::SimData* const Sim):
    Scheduler(Sim,"NeighbourListScheduler", NULL)
  { 
    dout << "Neighbour List Scheduler Algorithmn Loaded" << std::endl;
    operator<<(XML);
  }

  SNeighbourList::SNeighbourList(dynamo::SimData* const Sim, CSSorter* ns):
    Scheduler(Sim,"NeighbourListScheduler", ns)
  { dout << "Neighbour List Scheduler Algorithmn Loaded" << std::endl; }

  void 
  SNeighbourList::getParticleNeighbourhood(const Particle& part,
					   const nbHoodFunc& func) const
  {
#ifdef DYNAMO_DEBUG
    if (!std::tr1::dynamic_pointer_cast<GNeighbourList>(Sim->dynamics.getGlobals()[NBListID]))
      M_throw() << "Not a GNeighbourList!";
#endif

    //Grab a reference to the neighbour list
    const GNeighbourList& nblist(*static_cast<const GNeighbourList*>
				 (Sim->dynamics.getGlobals()[NBListID]
				  .get()));
  
    nblist.getParticleNeighbourhood(part, func);
  }
    
  void 
  SNeighbourList::getLocalNeighbourhood(const Particle& part, 
					const nbHoodFunc& func) const
  {
#ifdef DYNAMO_DEBUG
    if (!std::tr1::dynamic_pointer_cast<GNeighbourList>(Sim->dynamics.getGlobals()[NBListID]))
      M_throw() << "Not a GNeighbourList!";
#endif

    //Grab a reference to the neighbour list
    const GNeighbourList& nblist(*static_cast<const GNeighbourList*>
				 (Sim->dynamics.getGlobals()[NBListID]
				  .get()));

    //Add the local cell events
    nblist.getLocalNeighbourhood(part, func);
  }

  void 
  SNeighbourList::getParticleNeighbourhood(const Vector& vec, const nbHoodFunc2& func) const
  {
#ifdef DYNAMO_DEBUG
    if (!std::tr1::dynamic_pointer_cast<GNeighbourList>(Sim->dynamics.getGlobals()[NBListID]))
      M_throw() << "Not a GNeighbourList!";
#endif

    //Grab a reference to the neighbour list
    const GNeighbourList& nblist(*static_cast<const GNeighbourList*>
				 (Sim->dynamics.getGlobals()[NBListID]
				  .get()));
    
    nblist.getParticleNeighbourhood(vec, func);
  }

}
