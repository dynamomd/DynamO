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

#include <dynamo/dynamics/systems/nblistCompressionFix.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/dynamics/units/units.hpp>
#include <dynamo/dynamics/globals/neighbourList.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <dynamo/dynamics/NparticleEventData.hpp>
#include <dynamo/schedulers/scheduler.hpp>

#ifdef DYNAMO_DEBUG 
#include <boost/math/special_functions/fpclassify.hpp>
#endif

namespace dynamo {
  SysNBListCompressionFix::SysNBListCompressionFix(dynamo::SimData* nSim, double nGR, size_t nblistID):
    System(nSim),
    growthRate(nGR),
    cellID(nblistID)
  {
    sysName = "GlobalCellsCompressionHack";
    type = NON_EVENT;

    if (!std::tr1::dynamic_pointer_cast<GNeighbourList>(Sim->dynamics.getGlobals()[cellID]))
      M_throw() << "The ID passed to SysNBListCompressionFix isn't a GNeighbourList";
  }

  void
  SysNBListCompressionFix::initialise(size_t nID)
  {
    ID = nID;

    if (!std::tr1::dynamic_pointer_cast<GNeighbourList>(Sim->dynamics.getGlobals()[cellID]))
      M_throw() << "Have the globals been shuffled? The cellID is no longer a GNeighbourList.";
  
    GNeighbourList& nblist(dynamic_cast<GNeighbourList&>(*Sim->dynamics.getGlobals()[cellID]));

    initialSupportedRange = nblist.getMaxInteractionRange();
      
    dt = (nblist.getMaxSupportedInteractionLength() / initialSupportedRange - 1.0) 
      / growthRate - Sim->dSysTime;

    dout << "Compression Hack Loaded"
	 << "\nFor global " << nblist.getName()
	 << "\nCompression rate = " 
	 << growthRate / Sim->dynamics.units().unitTime()
	 << "\nSim Units compression rate = " << growthRate
	 << "\nMax length of interaction = " 
	 << nblist.getMaxSupportedInteractionLength() / Sim->dynamics.units().unitLength()
	 << "\nMaximum supported length = "
	 << nblist.getMaxSupportedInteractionLength() / Sim->dynamics.units().unitLength()
	 << "\nFirst halt scheduled for " 
	 << dt / Sim->dynamics.units().unitTime() << std::endl;
  }

  void
  SysNBListCompressionFix::runEvent() const
  {
    double locdt = dt;
  
#ifdef DYNAMO_DEBUG 
    if (boost::math::isnan(dt))
      M_throw() << "A NAN system event time has been found";
#endif
  
    Sim->dSysTime += locdt;
  
    Sim->ptrScheduler->stream(locdt);
  
    //dynamics must be updated first
    Sim->dynamics.stream(locdt);

    Sim->freestreamAcc += locdt;

    GNeighbourList& nblist(dynamic_cast<GNeighbourList&>
			   (*Sim->dynamics.getGlobals()[cellID]));
  
    dout << "Rebuilding the neighbour list named " << nblist.getName()
	 << "\nNColl = " << Sim->eventCount
	 << "\nSys t = " << Sim->dSysTime / Sim->dynamics.units().unitTime() << std::endl;
  
    nblist.setMaxInteractionRange(nblist.getMaxSupportedInteractionLength() * 1.1);
  
    dt = (nblist.getMaxSupportedInteractionLength()
	  / initialSupportedRange - 1.0) / growthRate - Sim->dSysTime;
  }

  void 
  SysNBListCompressionFix::fixNBlistForOutput()
  {
    GNeighbourList& nblist(dynamic_cast<GNeighbourList&>
			   (*Sim->dynamics.getGlobals()[cellID]));

    nblist.setMaxInteractionRange(initialSupportedRange * (1.0 + growthRate * Sim->dSysTime));
  }

}
