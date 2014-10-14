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

#include <dynamo/systems/nblistCompressionFix.hpp>

#include <dynamo/units/units.hpp>
#include <dynamo/globals/neighbourList.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>

namespace dynamo {
  SysNBListCompressionFix::SysNBListCompressionFix(dynamo::Simulation* nSim, double nGR, size_t nblistID):
    System(nSim),
    growthRate(nGR),
    cellID(nblistID)
  {
    sysName = "GlobalCellsCompressionHack";
    type = NON_EVENT;

    if (!std::dynamic_pointer_cast<GNeighbourList>(Sim->globals[cellID]))
      M_throw() << "The ID passed to SysNBListCompressionFix isn't a GNeighbourList";
  }

  void
  SysNBListCompressionFix::initialise(size_t nID)
  {
    ID = nID;

    if (!std::dynamic_pointer_cast<GNeighbourList>(Sim->globals[cellID]))
      M_throw() << "Have the globals been shuffled? The cellID is no longer a GNeighbourList.";
  
    GNeighbourList& nblist(dynamic_cast<GNeighbourList&>(*Sim->globals[cellID]));

    initialSupportedRange = nblist.getMaxInteractionRange();
      
    dt = (nblist.getMaxSupportedInteractionLength() / initialSupportedRange - 1.0) 
      / growthRate - Sim->systemTime;

    dout << "Compression Hack Loaded"
	 << "\nFor global " << nblist.getName()
	 << "\nCompression rate = " 
	 << growthRate / Sim->units.unitTime()
	 << "\nSim Units compression rate = " << growthRate
	 << "\nMax length of interaction = " 
	 << nblist.getMaxSupportedInteractionLength() / Sim->units.unitLength()
	 << "\nMaximum supported length = "
	 << nblist.getMaxSupportedInteractionLength() / Sim->units.unitLength()
	 << "\nFirst halt scheduled for " 
	 << dt / Sim->units.unitTime() << std::endl;
  }

  void
  SysNBListCompressionFix::runEvent()
  {
    Event event = getEvent();
  
#ifdef DYNAMO_DEBUG 
    if (std::isnan(event._dt))
      M_throw() << "A NAN system event time has been found";
#endif
  
    Sim->systemTime += event._dt;
  
    Sim->ptrScheduler->stream(event._dt);
  
    //dynamics must be updated first
    Sim->stream(event._dt);

    GNeighbourList& nblist(dynamic_cast<GNeighbourList&>(*Sim->globals[cellID]));
  
    dout << "Rebuilding the neighbour list named " << nblist.getName()
	 << "\nNColl = " << Sim->eventCount
	 << "\nSys t = " << Sim->systemTime / Sim->units.unitTime() << std::endl;
  
    nblist.setMaxInteractionRange(nblist.getMaxSupportedInteractionLength() * 1.1);
  
    dt = (nblist.getMaxSupportedInteractionLength()
	  / initialSupportedRange - 1.0) / growthRate - Sim->systemTime;

    NEventData SDat;

    Sim->_sigParticleUpdate(SDat);
    
    for (shared_ptr<OutputPlugin>& Ptr : Sim->outputPlugins)
      Ptr->eventUpdate(event, SDat); 
  }

  void 
  SysNBListCompressionFix::fixNBlistForOutput()
  {
    GNeighbourList& nblist(dynamic_cast<GNeighbourList&>
			   (*Sim->globals[cellID]));

    nblist.setMaxInteractionRange(initialSupportedRange * (1.0 + growthRate * Sim->systemTime));
  }

}
